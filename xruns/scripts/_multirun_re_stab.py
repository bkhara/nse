import argparse
import io
import os
import re
import shutil
import subprocess
import time
import fileinput

import libconf

from _exec_settings import *

egf = ExecGlobalFileNames
pas = ProgramArgumentStrings

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKCYANBOLD = '\033[1;96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    FAILBOLD = '\033[1;91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def get_args(parser):
    parser.add_argument(
        "--runtype", "-type",
        type=str,
        choices=[pas.PREP, pas.RUN, pas.DEV, pas.PROD, pas.PP, pas.CONS, pas.CHECKDIV],
        help=f"either {{ {pas.PREP} | {pas.RUN} | {pas.DEV} | {pas.PROD} | {pas.PP} | {pas.CONS} | {pas.CHECKDIV} }}",
        required=True,
    )
    return parser.parse_args()


def get_run_script_filename(arg):
    if arg.runtype == pas.DEV:
        return RunScript.cluster_dev
    if arg.runtype == pas.PROD:
        return RunScript.cluster_prod
    return RunScript.simple


def set_project_path_in_script(script_filename):
    """Keep the same ${PROJECT_DIR} replacement style as the old scripts."""
    if not os.path.exists(script_filename):
        return

    if not os.path.exists(PROJECT_PATH):
        raise FileNotFoundError(f"PROJECT_PATH = {PROJECT_PATH} does not exist.")

    backup_name_1 = f"{script_filename}.bak"
    backup_name_2 = script_filename.replace(".sh", "_backup.sh")
    backup_exists = os.path.exists(backup_name_1) or os.path.exists(backup_name_2)
    backup_ext = '' if backup_exists else '.bak'

    with fileinput.FileInput(script_filename, inplace=True, backup=backup_ext) as file:
        for line in file:
            print(line.replace("${PROJECT_DIR}", PROJECT_PATH), end='')

    if os.path.exists(backup_name_1):
        shutil.move(backup_name_1, backup_name_2)


def safe_case_float_string(value):
    """Use compact directory names, while remaining safe for decimal Re later."""
    return str(value).replace('.', 'p')


def case_dir_name(reynolds_number, stab_scheme):
    return f"re-{safe_case_float_string(reynolds_number)}-stab_{stab_scheme}"


def solve_time_filename(stab_scheme):
    return f"solve_times_stab_{stab_scheme}.txt"


def checkdiv_filename():
    return "divergence_summary.txt"


class ReStabRuns(object):
    """
    Generic Reynolds-number/stabilization sweep driver.

    This base class owns the common mechanics:
      - creating case directories
      - modifying/copying config.txt
      - copying supporting files
      - running/submitting jobs
      - consolidation/check-divergence utilities

    Case-specific subclasses should define:
      - rev
      - decide_time_span_and_step(re)
      - atomic_postprocess()
    """

    rev = []
    case_name = "generic"

    def __init__(self, arg):
        super(ReStabRuns, self).__init__()
        self.top_path = os.getcwd()
        print("in constructor, top path = ", self.top_path)

        self.config_template_path = os.path.join(self.top_path, "config.txt")
        with io.open(self.config_template_path) as f:
            self.config = libconf.load(f)

        self.arg = arg
        self.idata = AttrDict()
        self.idata.Re = None
        self.idata.T = None
        self.idata.dt = None
        self.idata.stab_scheme = None

        if not self.rev:
            raise ValueError(
                f"No Reynolds-number sweep is defined for {self.__class__.__name__}. "
                "Define a non-empty class variable `rev` in the case-specific subclass."
            )
        self.rev = list(self.rev)

        self.common_supporting_files = SupportingFiles.common
        self.specific_supporting_files = self.infer_specific_supporting_files()

        self.is_monitoring = (self.arg.runtype == pas.CHECKDIV)
        if not self.is_monitoring:
            with open(os.path.join(self.top_path, egf.ExecLogFile), "a") as f:
                f.write(f"--------------------------------\n")
                f.write(f"New {self.arg.runtype}\n")
                f.write(f"case = {self.case_name}\n")
                f.write(f"Re = {self.rev}\n")
                f.write(f"stab_scheme = {STAB_SCHEMES}\n")
                f.write(f"--------------------------------\n")

        if self.arg.runtype == pas.CONS:
            self.setup_consolidation_files()
        if self.arg.runtype == pas.CHECKDIV:
            with open(os.path.join(self.top_path, checkdiv_filename()), "w") as f:
                f.write("Re,stab_scheme,num_diverged\n")

    @staticmethod
    def decide_time_span_and_step(re_value):
        raise NotImplementedError(
            "Case-specific subclasses must implement decide_time_span_and_step(re_value)."
        )

    def infer_specific_supporting_files(self):
        files = []

        pcase = self.config.get('pcase', None)
        if pcase == FlowPastCylinder:
            files.extend(SupportingFiles.fpc)
        elif pcase == LidDrivenCavity:
            files.extend(SupportingFiles.ldc)

        # Also copy the actual mesh named by this newer config.
        mesh_config = self.config.get('mesh_config', {})
        mesh_file = mesh_config.get('mesh_file', None)
        if mesh_file is not None:
            files.append(mesh_file)

        # Preserve order and remove duplicates.
        return list(dict.fromkeys(files))

    def copy_supporting_files(self):
        for filename in self.common_supporting_files:
            src = os.path.join(self.top_path, filename)
            if os.path.exists(src):
                shutil.copy(src, ".")

        for filename in self.specific_supporting_files:
            src = os.path.join(self.top_path, filename)
            if os.path.exists(src):
                shutil.copy(src, ".")

    def modify_config(self):
        self.config['Re'] = float(self.idata.Re)

        if 'time_marching' not in self.config:
            self.config['time_marching'] = {}

        self.config['time_marching']['t_max'] = float(self.idata.T)
        self.config['time_marching']['dt'] = float(self.idata.dt)

        if 'method_config' not in self.config:
            self.config['method_config'] = {}
        self.config['method_config']['stab_scheme'] = self.idata.stab_scheme

        # Keep this explicit so the config is internally consistent when toggling
        # between no stabilization and SUPG/PSPG.
        if 'sups_config' not in self.config:
            self.config['sups_config'] = {}
        use_sups = (self.idata.stab_scheme == "sups")
        self.config['sups_config']['use_supg'] = use_sups
        self.config['sups_config']['use_pspg'] = USE_PSPG

        nsteps = int(round(float(self.idata.T) / float(self.idata.dt)))
        self.config['file_write_freq'] = max(1, int(round(nsteps / FIXED_NUM_FILES_TO_WRITE)))

    def copy_config(self):
        with io.open('./config.txt', 'w') as f:
            libconf.dump(self.config, f)

    def create_case_dir(self, dir_string):
        case_dir = os.path.join(".", dir_string)
        if not os.path.exists(case_dir):
            os.makedirs(case_dir)
        return case_dir

    def atomic_preprocess(self):
        self.modify_config()
        self.copy_config()
        self.copy_supporting_files()

    def atomic_run(self):
        run_script = RunScript.simple
        if not os.path.exists(run_script):
            raise FileNotFoundError(f"{run_script} not found")
        subprocess.run(["bash", run_script])

    def atomic_prod_run(self):
        job_script = RunScript.cluster_prod
        if not os.path.exists(job_script):
            raise FileNotFoundError(f"{job_script} not found")
        result = subprocess.run(["sbatch", job_script], capture_output=True, text=True)
        jobId = re.findall(r'Submitted batch job (\d+)', result.stdout)
        with open("submitted_job.txt", "a") as file:
            file.write(jobId[0] + "\n" if jobId else result.stdout + result.stderr + "\n")

    def atomic_dev_run(self):
        job_script = RunScript.cluster_dev
        if not os.path.exists(job_script):
            raise FileNotFoundError(f"{job_script} not found")
        result = subprocess.run(["sbatch", job_script], capture_output=True, text=True)
        jobId = re.findall(r'Submitted batch job (\d+)', result.stdout)
        with open("submitted_job.txt", "a") as file:
            file.write(jobId[0] + "\n" if jobId else result.stdout + result.stderr + "\n")

    def atomic_postprocess(self):
        print(
            f"{bcolors.WARNING}No postprocessing implemented for "
            f"{self.__class__.__name__}: "
            f"Re={self.idata.Re}, stab_scheme={self.idata.stab_scheme}{bcolors.ENDC}"
        )

    def setup_consolidation_files(self):
        for stab_scheme in STAB_SCHEMES:
            with open(os.path.join(self.top_path, solve_time_filename(stab_scheme)), "w") as file:
                file.write("Re,solve_time\n")

    def atomic_consolidation(self):
        Re = self.idata.Re
        stab_scheme = self.idata.stab_scheme
        filename = solve_time_filename(stab_scheme)

        pattern = r'\[TIME\] Solve-time \(global_total_sec avg\):\s+([0-9.eE+-]+)'
        out_file_path = egf.Outfile
        if os.path.isfile(out_file_path):
            with open(out_file_path, 'r') as outfile:
                for line in outfile:
                    match = re.search(pattern, line)
                    if match:
                        with open(os.path.join(self.top_path, filename), 'a') as file:
                            file.write(f"{Re},{match.group(1)}\n")
                        break

    def atomic_check_divergent_behavior(self):
        Re = self.idata.Re
        stab_scheme = self.idata.stab_scheme

        pattern = r'DIVERGED'
        out_file_path = egf.Outfile
        num_diverged = 0
        if os.path.isfile(out_file_path):
            with open(out_file_path, 'r') as outfile:
                content = outfile.read()
            num_diverged = len(re.findall(pattern, content))

        if num_diverged > 0:
            print(f"{bcolors.FAIL}Re={Re}, stab_scheme={stab_scheme}: DIVERGED x {num_diverged}{bcolors.ENDC}")
        else:
            print(f"{bcolors.OKGREEN}Re={Re}, stab_scheme={stab_scheme}: no DIVERGED pattern found{bcolors.ENDC}")

        with open(os.path.join(self.top_path, checkdiv_filename()), "a") as f:
            f.write(f"{Re},{stab_scheme},{num_diverged}\n")

    def action(self, case_string):
        cwd = os.getcwd()
        case_dir = self.create_case_dir(case_string)
        os.chdir(case_dir)

        print(f"{bcolors.WARNING}{case_string}{bcolors.ENDC}: {os.getcwd()}")
        if not self.is_monitoring:
            with open(os.path.join(self.top_path, egf.ExecLogFile), "a") as f:
                f.write(f"{case_string}\n")

        if self.arg.runtype == pas.PREP:
            self.atomic_preprocess()
        elif self.arg.runtype == pas.RUN:
            self.atomic_preprocess()
            self.atomic_run()
        elif self.arg.runtype == pas.DEV:
            self.atomic_preprocess()
            self.atomic_dev_run()
        elif self.arg.runtype == pas.PROD:
            self.atomic_preprocess()
            self.atomic_prod_run()
        elif self.arg.runtype == pas.PP:
            self.atomic_postprocess()
        elif self.arg.runtype == pas.CONS:
            self.atomic_consolidation()
        elif self.arg.runtype == pas.CHECKDIV:
            self.atomic_check_divergent_behavior()

        os.chdir(cwd)

    def run_all_cases(self):
        for re_value in self.rev:
            self.idata.Re = re_value

            T, dt = self.decide_time_span_and_step(re_value)
            self.idata.T = T
            self.idata.dt = dt

            for stab_scheme in STAB_SCHEMES:
                self.idata.stab_scheme = stab_scheme
                self.action(case_dir_name(re_value, stab_scheme))


def run_runner(runner_class):
    parser = argparse.ArgumentParser(
        prog=f'{runner_class.__name__}',
        usage='python %(prog)s [options]'
    )
    arg = get_args(parser)

    start_time = time.perf_counter()

    run_script_filename = get_run_script_filename(arg)
    set_project_path_in_script(run_script_filename)

    runner = runner_class(arg)
    runner.run_all_cases()

    end_time = time.perf_counter()
    print(
        bcolors.FAILBOLD,
        f"Elapsed time = {(end_time-start_time):.2f} sec ({(end_time-start_time) / 60:.2f} min)",
        bcolors.ENDC,
    )


if __name__ == "__main__":
    # The generic base class is intentionally not case-aware anymore.
    # Use multirun_cases.py as the normal entry point.
    print(
        f"{bcolors.WARNING}Use `python multirun_cases.py -type <prep|run|dev|prod|pp|cons|checkdiv>` "
        f"so the correct case-specific subclass can be selected from config.txt.{bcolors.ENDC}"
    )
