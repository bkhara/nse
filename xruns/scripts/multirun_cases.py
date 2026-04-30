import argparse
import glob
import io
import os
import subprocess
import time

import libconf

from _exec_settings import *
from _multirun_re_stab import (
    ReStabRuns,
    bcolors,
    get_args,
    get_run_script_filename,
    set_project_path_in_script,
)


class LDCReStabRuns(ReStabRuns):
    case_name = "ldc"

    rev = [
        100,
        400,
        1000,
        3200,
        5000,
        7500,
        10000,
    ]

    @staticmethod
    def decide_time_span_and_step(re_value):
        totalT = 50
        dt = 1

        if re_value < 5000:
            totalT = 5000  # 250 is sufficient
            dt = 1.0
        elif 5000 <= re_value < 7000:
            totalT = 5000  # 1500 is sufficient
            dt = 0.5
        elif 7000 <= re_value:
            totalT = 5000
            dt = 0.1

        return totalT, dt

    def atomic_postprocess(self):
        Re = self.idata.Re
        stab_scheme = self.idata.stab_scheme

        print(
            f"{bcolors.OKCYAN}LDC postprocess: "
            f"Re={Re}, stab_scheme={stab_scheme}{bcolors.ENDC}"
        )

        pvtu_file = self.find_latest_pvtu()
        if pvtu_file is None:
            print(
                f"{bcolors.WARNING}No Cycle*/data.pvtu found for "
                f"Re={Re}, stab_scheme={stab_scheme}{bcolors.ENDC}"
            )
            return

        self.extract_ldc_centerlines(pvtu_file)

    @staticmethod
    def find_latest_pvtu():
        candidates = glob.glob("Cycle*/data.pvtu")
        if not candidates:
            return None

        def cycle_number(path):
            dirname = os.path.basename(os.path.dirname(path))
            try:
                return int(dirname.replace("Cycle", ""))
            except ValueError:
                return -1

        return max(candidates, key=cycle_number)

    def extract_ldc_centerlines(self, pvtu_file):
        import numpy as np
        import pandas as pd
        import pyvista as pv

        Re = self.idata.Re
        stab_scheme = self.idata.stab_scheme

        mesh = pv.read(pvtu_file)
        vel_name = self.find_velocity_array_name(mesh)

        npts = 300

        def sample_line(points, outname):
            sampled = pv.PolyData(points).sample(mesh)
            u = sampled[vel_name]

            df = pd.DataFrame({
                "Re": Re,
                "x": points[:, 0],
                "y": points[:, 1],
                "z": points[:, 2],
                "ux": u[:, 0],
                "uy": u[:, 1],
            })

            df.to_csv(outname, index=False)
            print(f"{bcolors.OKGREEN}Wrote {outname}{bcolors.ENDC}")

        # Vertical centerline: x = 0.5, y in [0, 1].
        y = np.linspace(0.0, 1.0, npts)
        vertical_points = np.column_stack([
            0.5 * np.ones_like(y),
            y,
            np.zeros_like(y),
            ])

        vertical_out = f"re-{Re}-stab_{stab_scheme}_vertical_centerline.csv"
        sample_line(vertical_points, vertical_out)

        # Horizontal centerline: y = 0.5, x in [0, 1].
        x = np.linspace(0.0, 1.0, npts)
        horizontal_points = np.column_stack([
            x,
            0.5 * np.ones_like(x),
            np.zeros_like(x),
            ])

        horizontal_out = f"re-{Re}-stab_{stab_scheme}_horizontal_centerline.csv"
        sample_line(horizontal_points, horizontal_out)

    @staticmethod
    def find_velocity_array_name(mesh):
        candidates = [
            "velocity",
            "Velocity",
            "u",
            "U",
            "vel",
        ]

        for name in candidates:
            if name in mesh.point_data:
                return name

        vector_arrays = []
        for name, arr in mesh.point_data.items():
            if len(arr.shape) == 2 and arr.shape[1] >= 2:
                vector_arrays.append(name)

        if len(vector_arrays) == 1:
            return vector_arrays[0]

        raise RuntimeError(
            "Could not infer velocity array name. "
            f"Available point arrays are: {list(mesh.point_data.keys())}"
        )


class FPCReStabRuns(ReStabRuns):
    case_name = "fpc"

    rev = [
        1,
        15,
        30,
        60,
        100,
        150,
        200,
        250,
        300,
        400,
    ]

    @staticmethod
    def decide_time_span_and_step(re_value):
        if re_value <= 50:
            totalT = 100
            dt = 0.2
        elif 50 < re_value <= 100:
            totalT = 150
            dt = 0.1
        elif 100 < re_value <= 200:
            totalT = 200
            dt = 0.05
        elif 200 < re_value <= 300:
            totalT = 250
            dt = 0.03
        elif 300 < re_value:
            totalT = 300
            dt = 0.01
        else:
            totalT = 300
            dt = 0.01

        return totalT, dt

    def atomic_postprocess(self):
        Re = self.idata.Re
        stab_scheme = self.idata.stab_scheme

        print(
            f"{bcolors.OKCYAN}FPC postprocess: "
            f"Re={Re}, stab_scheme={stab_scheme}{bcolors.ENDC}"
        )

        # Keep this compatible with your existing FPC helper.
        if os.path.exists("analyse_forces.py"):
            subprocess.run(["python", "analyse_forces.py"], check=False)
        else:
            print(
                f"{bcolors.WARNING}analyse_forces.py not found for "
                f"Re={Re}, stab_scheme={stab_scheme}{bcolors.ENDC}"
            )


def make_runner(arg):
    with io.open("config.txt") as f:
        config = libconf.load(f)

    pcase = config.get("pcase", None)

    if pcase == LidDrivenCavity:
        return LDCReStabRuns(arg)

    if pcase == FlowPastCylinder:
        return FPCReStabRuns(arg)

    raise ValueError(
        f"No case-specific multirun subclass is available for pcase={pcase}."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="multirun_cases.py",
        usage="python %(prog)s [options]",
    )
    arg = get_args(parser)

    start_time = time.perf_counter()

    run_script_filename = get_run_script_filename(arg)
    set_project_path_in_script(run_script_filename)

    runner = make_runner(arg)
    runner.run_all_cases()

    end_time = time.perf_counter()
    print(
        bcolors.FAILBOLD,
        f"Elapsed time = {(end_time-start_time):.2f} sec ({(end_time-start_time) / 60:.2f} min)",
        bcolors.ENDC,
    )
