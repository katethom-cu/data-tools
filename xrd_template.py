#!/usr/bin/env python3
import matplotlib.pyplot as plt
from xrd import XRDSample, xrd_multiplot

plt.rcParams["text.usetex"] = True


def main():
    samples = {
        "CR-001-A": XRDSample.from_files(
            "./Cr-001-post_I4_20250516.ASC", "./cr-001-post-data.xlsx"
        ).normalize(),
        "CR-001-B": XRDSample.from_files(
            "./C1-001A_I14_20250523.ASC", "./Cr-001A-data.xlsx"
        ).normalize(),
        "CR-001-C": XRDSample.from_files(
            "./C1-001B_I15_20250523.ASC", "./Cr-001-B-data.xlsx"
        ).normalize(),
    }

    xrd_multiplot(samples)

    plt.show()


if __name__ == "__main__":
    main()
