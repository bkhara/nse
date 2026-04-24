#!/usr/bin/env python3
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


def plot_integer_hist(counts, xlabel, title, outname, logy=False):
    counts = counts.astype(int)

    vmin = counts.min()
    vmax = counts.max()
    bins = range(vmin, vmax + 2)

    plt.figure(figsize=(7, 5))
    plt.hist(counts, bins=bins, align="left", rwidth=0.85)
    plt.xticks(range(vmin, vmax + 1))
    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.title(title)

    if logy:
        plt.yscale("log")

    plt.tight_layout()
    plt.savefig(outname, dpi=200)
    plt.close()


def format_top_n(series, n=5):
    vals = series.sort_values(ascending=False).head(n).tolist()
    return ", ".join(str(int(v)) for v in vals)


def main():
    stag_file = Path("stag_history.txt")
    pg_file = Path("pg_history.txt")
    out_file = Path("avgiters.txt")

    if not stag_file.is_file():
        print(f"Error: required file not found: {stag_file}")
        return 1

    lines = []

    # ------------------------------------------
    # Read staggered history
    # ------------------------------------------
    df_stag = pd.read_csv(stag_file)

    if not {"t", "stag_it"}.issubset(df_stag.columns):
        print("Error: stag_history.txt must contain at least columns: t, stag_it")
        return 1

    stem = "history"

    # -------------------------------------------------
    # Staggered iteration counts per time step
    # count = max(stag_it) + 1 for each t
    # -------------------------------------------------
    stag_counts = (df_stag.groupby("t")["stag_it"].max() + 1).astype(int)

    n_time_steps = len(stag_counts)

    total_stag = stag_counts.sum()
    avg_stag = stag_counts.mean()
    max_stag = stag_counts.max()
    min_stag = stag_counts.min()
    t_of_max_stag = stag_counts.idxmax()
    t_of_min_stag = stag_counts.idxmin()
    top5_stag = format_top_n(stag_counts, n=5)

    lines.append(f"Number of time steps: {n_time_steps}")
    lines.append(f"Total staggered iterations: {total_stag}")
    lines.append(f"Average staggered iterations per time step = {avg_stag:.6f}")
    lines.append(f"Maximum staggered iterations in a time step = {max_stag}")
    lines.append(f"Minimum staggered iterations in a time step = {min_stag}")
    lines.append(f"Top 5 staggered-iteration counts over time steps = {top5_stag}")
    lines.append(f"Time step where maximum staggered iteration happens = {t_of_max_stag:.16g}")
    lines.append(f"Time step where minimum staggered iteration happens = {t_of_min_stag:.16g}")

    stag_hist_file = f"{stem}_stag_hist.png"
    plot_integer_hist(
        stag_counts,
        xlabel="Number of staggered iterations per time step",
        title="Histogram of staggered iterations",
        outname=stag_hist_file,
        logy=False,
    )
    lines.append(f"Saved staggered histogram to: {stag_hist_file}")

    # ------------------------------------------
    # Read PG history only if available
    # ------------------------------------------
    if pg_file.is_file():
        df_pg = pd.read_csv(pg_file)

        if not {"t", "stag_it", "k"}.issubset(df_pg.columns):
            lines.append("Warning: pg_history.txt exists, but does not contain columns: t, stag_it, k")
        else:
            # -------------------------------------------------
            # PG iteration counts per staggered step
            # count = max(k) + 1 for each (t, stag_it)
            # -------------------------------------------------
            pg_counts = (df_pg.groupby(["t", "stag_it"])["k"].max() + 1).astype(int)

            avg_pg = pg_counts.mean()
            max_pg = pg_counts.max()
            min_pg = pg_counts.min()
            t_stag_of_max_pg = pg_counts.idxmax()   # tuple: (t, stag_it)
            t_stag_of_min_pg = pg_counts.idxmin()   # tuple: (t, stag_it)
            top5_pg = format_top_n(pg_counts, n=5)

            lines.append("")
            lines.append(f"Average PG iterations per staggered iteration = {avg_pg:.6f}")
            lines.append(f"Maximum PG iterations in a staggered step = {max_pg}")
            lines.append(f"Minimum PG iterations in a staggered step = {min_pg}")
            lines.append(f"Top 5 PG-iteration counts over staggered steps = {top5_pg}")
            lines.append(
                "Location where maximum PG iteration happens = "
                f"(t={t_stag_of_max_pg[0]:.16g}, stag_it={t_stag_of_max_pg[1]})"
            )
            lines.append(
                "Location where minimum PG iteration happens = "
                f"(t={t_stag_of_min_pg[0]:.16g}, stag_it={t_stag_of_min_pg[1]})"
            )

            pg_hist_file = f"{stem}_pg_hist.png"
            plot_integer_hist(
                pg_counts,
                xlabel="Number of PG iterations per staggered step",
                title="Histogram of PG iterations",
                outname=pg_hist_file,
                logy=True,
            )
            lines.append(f"Saved PG histogram to: {pg_hist_file}")
    else:
        lines.append("")
        lines.append("pg_history.txt not found. Skipping PG iteration statistics and histogram.")

    out_file.write_text("\n".join(lines) + "\n")
    print(f"Wrote output to: {out_file}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())