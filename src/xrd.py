import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from style import MARKERS, LINE_COLOURS

plt.rc("text", usetex=True)
plt.rc("text.latex", preamble=r"\usepackage{mhchem}")



def read_ascii(filename):
    with open(filename, "r") as f:
        lines = np.array([line.split(" ") for line in f], dtype=float)
        return np.hsplit(lines, 3)


class Compound:
    def __init__(self, name: str, peak_theta2: np.ndarray, peak_intensity: np.ndarray):
        self.name = name
        self.peak_theta2 = peak_theta2
        self.peak_intensity = peak_intensity

    def normalize(self, max=None):
        if max is None:
            max = np.max(self.peak_intensity)

        norm_intensities = self.peak_intensity / max

        return Compound(self.name, self.peak_theta2, norm_intensities)

    def offset(self, amount):
        return Compound(self.name, self.peak_theta2, self.peak_intensity + amount)

    def label_peaks(self, ax, marker="v", peak_label_offset=0.03):
        max = np.max(self.peak_intensity)
        return ax.scatter(
            self.peak_theta2,
            self.peak_intensity + max * peak_label_offset,
            label=f"\\ce{{{self.name.replace(" ", "")}}}",
            marker=marker,
            color="black",
        )


class XRDSample:
    def __init__(
        self, theta2: np.ndarray, intensity: np.ndarray, compounds: list[Compound]
    ):
        self.theta2 = theta2
        self.intensity = intensity
        self.compounds = compounds

    @classmethod
    def from_files(cls, ascii_filename, excel_filename):
        theta2, intensity, _ = read_ascii(ascii_filename)
        meta = pd.read_excel(excel_filename, sheet_name=None)

        compounds = []

        for _, row in meta["info"][["Chemical Formula", "Ref. Code"]].iterrows():
            name = row["Chemical Formula"]
            ref = row["Ref. Code"]

            sample_peaks = meta["sample"]
            matched_peak_positions = np.array(
                sample_peaks[sample_peaks["Matched by"].str.contains(ref, na=False)]["Pos. [°2Th.]"]
            )
            indices = [
                np.abs(theta2 - peak).argmin() for peak in matched_peak_positions
            ]

            # fuck this shit!!!
            # why the hell isn't the peak angle at the maxima. i hate this.
            fuck_this_shit = 10
            peak_intensities = np.array(
                [
                    np.max(intensity[peak - fuck_this_shit : peak + fuck_this_shit])
                    for peak in indices
                ]
            )

            compounds.append(
                Compound(str(name), matched_peak_positions, peak_intensities)
            )

        return cls(theta2, intensity, compounds)

    def normalize(self):
        max = np.max(self.intensity)
        normalized_intensity = self.intensity / max
        normalized_compounds = [compound.normalize(max) for compound in self.compounds]

        return XRDSample(self.theta2, normalized_intensity, normalized_compounds)

    def offset(self, amount):
        offset_intensity = self.intensity + amount

        offset_compounds = [compound.offset(amount) for compound in self.compounds]

        return XRDSample(self.theta2, offset_intensity, offset_compounds)

    def plot(self, title=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        fig.tight_layout()
        ax.set_ylabel("Intensity (A.U.)")
        ax.set_yticks([])
        ax.set_xlabel("2$\\theta$ (degrees)")
        line = ax.plot(self.theta2, self.intensity, label=title)

        scatters = [
            compound.label_peaks(ax, marker=MARKERS[i])
            for i, compound in enumerate(self.compounds)
        ]

        ax.legend(handles=scatters, loc="lower right")

        return line

    def max_intensity(self):
        return np.max(self.intensity)


def xrd_multiplot(
    samples: dict[str, XRDSample],
    title: str | None = None,
    margin: float = 0.1,
    peak_label_offset: float = 0.01,
    peak_label_gap: float = 0.04
):
    spacing = 1 + margin

    _, ax = plt.subplots(figsize=(12 * 0.6, 8 * 0.6))
    if title is not None:
        ax.set_title(title)
    ax.set_ylabel("Intensity (A.U.)")
    ax.set_yticks([])
    ax.set_xlabel("2$\\theta$ (degrees)")

    max_intensity = np.max([sample.max_intensity() for sample in samples.values()])

    abs_peak_offset = peak_label_offset * max_intensity
    margin = peak_label_gap * max_intensity

    # plotting lines

    offset_samples = {
        name: sample.offset((len(samples) - i - 1) * spacing * max_intensity)
        for i, (name, sample) in enumerate(samples.items())
    }

    lines = [
        ax.plot(
            sample.theta2,
            sample.intensity,
            label=(name),
            zorder=1,
            color=LINE_COLOURS[i],
        )
        for i, (name, sample) in enumerate(offset_samples.items())
    ]

    ax.add_artist(
        ax.legend(
            handles=[item for sublist in lines for item in sublist], loc="upper left"
        )
    )

    # collating compounds
    # hate this code
    compounds = {}
    for sample in offset_samples.values():
        for compound in sample.compounds:
            if compound.name not in compounds.keys():
                compounds[compound.name] = {"theta2": [], "intensity": []}
            compounds[compound.name]["theta2"] = np.append(
                compounds[compound.name]["theta2"], compound.peak_theta2
            )
            compounds[compound.name]["intensity"] = np.append(
                compounds[compound.name]["intensity"], compound.peak_intensity
            )

    # plotting compounds
    # there's probably a pretty functional way of doing this. I do not know it
    scatter = []
    ive_lost_the_plot = []
    for i, (name, value) in enumerate(compounds.items()):
        for angle in value["theta2"]:
            ive_lost_the_plot.append(angle)

        intensities_offest = [
            intensity + margin * ive_lost_the_plot.count(angle)
            for angle, intensity in zip(value["theta2"], value["intensity"])]
        scatter.append(ax.scatter(
            value["theta2"],
            intensities_offest + abs_peak_offset,
            label=f"\\ce{{{name.replace(" ", "")}}}",
            color="black",
            marker=MARKERS[i],
            zorder=2,
        ))

    ax.legend(handles=scatter, loc="upper right")
