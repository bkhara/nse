import os
import shutil
import argparse
from pathlib import Path

RESULT_DIR = "./results"
CASES_HOME = "./cases"
SCRIPTS_HOME = "./scripts"

class BColors:
    red = '\033[91m'
    green = '\033[92m'
    yellow = '\033[93m'
    blue = '\033[94m'
    magenta = '\033[95m'
    cyan = '\033[96m'
    white = '\033[97m'
    Endc = '\033[0m'
    Bold = '\033[1m'
    Underline = '\033[4m'
cl = BColors

class CaseStrings:
    mms2 = "mms"
    fpc2 = "fpc"
    ldc2 = "ldc"
cs = CaseStrings

possible_scripts_with_project_dir_variable = [
    "_exec_settings.py",
    "runex.sh",
    "job.sh",
    "nova.sh",
    "oscar.sh"
]

multirun_scripts = [
    "_exec_settings.py",
    "_multirun_re_stab.py",
    "multirun_cases.py"
]

def determine_destination_path(dest_dir):
    counter = 1
    dest_dir_mod = dest_dir
    while os.path.exists(dest_dir_mod):
        print(f"{cl.Bold}{cl.magenta}Warning:: path already exists: {dest_dir_mod}", cl.Endc)
        dest_dir_mod = f"{dest_dir}-{counter}"
        counter += 1
    print(f"{cl.Bold}{cl.cyan}Setup initiated at: {dest_dir_mod}{cl.Endc}")
    return dest_dir_mod

def determine_project_path():
    # this fine (i.e., initiate_run.py) is supposed to be in the project/xruns directory
    # so, deducing the project path from the current path
    project_path = os.path.dirname(os.path.normpath(os.getcwd()))
    return project_path

def update_project_path_in_runex(project_path, runex_path):
    p = Path(runex_path)
    p.write_text(p.read_text().replace("${PROJECT_DIR}", project_path))

class CopyScripts:
    def copy_restart_file(self, arg, destination_dir : str):
        restart_script = os.path.join(SCRIPTS_HOME, "restart.py")
        shutil.copy(restart_script, destination_dir)

    def check_multirun_eligibility(self, arg):
        assert_flag = arg.casename == cs.fpc2 or arg.casename == cs.ldc2
        assert assert_flag, f"{cl.red}{arg.casename} does not have multirun script{cl.Endc}"

    def copy_multirun_files(self, arg, destination_dir : str):
        self.check_multirun_eligibility(arg)
        for file in multirun_scripts:
            srcfullpath = os.path.join(SCRIPTS_HOME, file)
            shutil.copy(srcfullpath, destination_dir)
        # multirun_file = os.path.join(SCRIPTS_HOME, "multirun_re_stab.py")
        # exec_file = os.path.join(SCRIPTS_HOME, "exec_settings.py")
        # shutil.copy(multirun_file, destination_dir)
        # shutil.copy(exec_file, destination_dir)

    def copy_scripts(self, arg, destination_dir : str):
        if arg.restart:
            self.copy_restart_file(arg, destination_dir)
        if arg.multi:
            self.copy_multirun_files(arg, destination_dir)

def get_args(parser):
    parser.add_argument("casename", type=str, help=f"Should be one of [ {cs.mms2} | {cs.fpc2} | {cs.ldc2} ] ]")
    parser.add_argument("-s", "--suffix", type=str, help=f"A suffix string that specifies the particular run", default = "")
    parser.add_argument("-e", "--exec", help=f"Copy the exec run scripts", action='store_true')
    parser.add_argument("-m", "--multi", help=f"Copy the multi run scripts", action='store_true')
    parser.add_argument("-r", "--restart", help=f"Copy the restart script", action='store_true')
    parser.add_argument("-nco", "--no_create_out", help=f"Do not create an empty out.txt file", action='store_true')
    parser.add_argument("-ppath", "--project_path", help=f"Current path of the project directory", default="")
    parser.add_argument("-dest", "--dest_path", help=f"Current path of the project directory", default="")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='initiate_run', usage=f"{cl.Bold}{cl.red}python %(prog)s <casename> [options]{cl.Endc}")
    arg = get_args(parser)
    print(arg)

    # PROJECT AND CASE SOURCE PATHS
    if not arg.project_path:
        # determine project path
        project_path = determine_project_path()
        print(f"Project path: {project_path} (deduced)")
    else:
        project_path = arg.project_path
        print(f"Project path: {project_path} (user input)")
    case_source_dir = os.path.join(CASES_HOME, arg.casename)
    print(case_source_dir)
    if not os.path.exists(case_source_dir):
        raise ValueError(f"{cl.Bold}{cl.red}Error: Invalid case: {case_source_dir}{cl.Endc}")

    # DESTINATION PATH
    if not arg.dest_path:
        DEST_DIR_TOP = os.path.join(RESULT_DIR, arg.casename)
        print(f"Dest location: {DEST_DIR_TOP} (default)")
    else:
        DEST_DIR_TOP = arg.dest_path
        print(f"Dest location: {DEST_DIR_TOP} (user input)")
    if not os.path.exists(DEST_DIR_TOP):
        os.makedirs(DEST_DIR_TOP)
    destination_dir = os.path.join(DEST_DIR_TOP, arg.suffix)
    destination_dir = determine_destination_path(destination_dir)

    # create the object for copying scripts
    cps = CopyScripts()
    # check if the right parameters have been passed for multirun
    if arg.multi:
        cps.check_multirun_eligibility(arg)
    # finally copy the directory structure and other scripts if applicable
    shutil.copytree(case_source_dir, destination_dir)
    cps.copy_scripts(arg, destination_dir)

    # update the path of the project in the run scripts
    for filename in possible_scripts_with_project_dir_variable:
        full_file_path = os.path.join(destination_dir, filename)
        if os.path.exists(full_file_path):
            update_project_path_in_runex(project_path, full_file_path)

    # create an empty "out.txt" file in the destination_dir
    if not arg.no_create_out:
        outfile = "out.txt"
        with open(os.path.join(destination_dir, outfile), 'w') as f:
            pass
