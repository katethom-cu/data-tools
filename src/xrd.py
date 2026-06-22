from typing import Tuple
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from style import MARKERS, LINE_COLOURS
import re
from enum import Enum, auto
from uncertainties import ufloat, UFloat

import scienceplots
plt.style.use(["science"])
plt.rc("text", usetex=True)
plt.rc("text.latex", preamble=r"\usepackage{mhchem}\usepackage{siunitx}")

class PhaseLabel(Enum):
    NONE = auto()
    LINES = auto()
    MARKERS = auto()

def read_ascii(filename):
    with open(filename, "r") as f:
        lines = np.array([line.split(" ") for line in f], dtype=float)
        return lines.T

class Compound:
    def __init__(self, name: str, peak_theta2: np.ndarray, peak_intensity: np.ndarray, qty: UFloat = 0.0):
        # sorted by significance. makes filtering easier later
        pair = zip(peak_theta2, peak_intensity)
        by_significance = np.array(sorted(pair, key=lambda x: x[1], reverse=True))
        self.name = name
        self.peak_theta2 = by_significance[:, 0]
        self.peak_intensity = by_significance[:, 1]
        self.qty = qty


    def normalize(self, max=None):
        if max is None:
            max = np.max(self.peak_intensity)

        norm_intensities = self.peak_intensity / max

        return Compound(self.name, self.peak_theta2, norm_intensities, self.qty)

    def subtract_bg(self, bg_theta2, bg_intensity):
        indices = [
            np.abs(bg_theta2 - peak).argmin() for peak in self.peak_theta2
        ]
        new_intensities = [
            inten - bg_intensity[i]
            for i, inten in zip(indices, self.peak_intensity)
        ]
        return Compound(
            self.name, self.peak_theta2, np.array(new_intensities), self.qty)

    def offset(self, amount):
        return Compound(self.name, self.peak_theta2, self.peak_intensity + amount)

    def label_peaks(self, ax, marker="v", peak_label_offset=0.04):
        max = np.max(self.peak_intensity)
        return ax.scatter(
            self.peak_theta2,
            self.peak_intensity + max * peak_label_offset,
            label=f"\\ce{{{self.name.replace(" ", "")}}}",
            marker=marker,
            color="black",
        )

    def __str__(self) -> str:
        return f"Name: {self.name}, {self.qty}%\n{pd.DataFrame({"2 Theta": self.peak_theta2, "Intensity": self.peak_intensity})}"

    def __repr__(self) -> str:
        return str(self)


class XRDSample:
    def __init__(
        self, theta2: np.ndarray, intensity: np.ndarray, compounds: list[Compound]
    ):
        self.theta2 = theta2
        self.intensity = intensity
        self.compounds = compounds

    @classmethod
    def from_profex(cls, raw_data, peaks, composition=None, globals=None):
        angle, intensity = read_ascii(raw_data)

        peaks = pd.read_csv(peaks,
                            dtype={
                                "Phase": str,
                                "h": np.int32,
                                "k": np.int32,
                                "l": np.int32,
                                "Angle (°2?)": np.float64,
                                "d (nm)": np.float64,
                                "Intensity (deg*cts)": np.float64,
                                "Rel. intensity (%)": np.float64,
                                "B1 (1/nm)": np.float64,
                                "B2 (1/nm)": np.float64,
                            })

        compound_names = peaks["Phase"].unique()

        # no idea what the intensity in that csv is supposed to represent, so
        # instead i'm just going to use the angle to look up the intensity at
        # the closest point in the raw data.
        # unreadable! but hey it's technically a one liner so that means its good right?
        if globals != None:
            globals = pd.read_csv(globals, dtype={
                    "Parameter, Goal": str,
                    "Value": np.float64,
                    "ESD": np.float64
                })
            compounds = [
                Compound(compound_name,
                    peaks.loc[peaks["Phase"] == compound_name]["Angle (°2?)"],
                    [intensity[(np.abs(angle - peak_angle)).argmin()]
                         for peak_angle in
                         peaks.loc[peaks["Phase"] == compound_name]["Angle (°2?)"]],
                    100 * ufloat(globals.loc[globals["Parameter, Goal"] == f"Q{compound_name}"]["Value"].iloc[0],
                                 globals.loc[globals["Parameter, Goal"] == f"Q{compound_name}"]["ESD"].iloc[0]))
                for compound_name in compound_names]

            print(compounds)

        elif composition == None:
            compounds = [
                Compound(compound_name,
                    peaks.loc[peaks["Phase"] == compound_name]["Angle (°2?)"],
                    [intensity[(np.abs(angle - peak_angle)).argmin()]
                         for peak_angle in
                         peaks.loc[peaks["Phase"] == compound_name]["Angle (°2?)"]])
                for compound_name in compound_names]
        else:
            qtys = pd.read_csv(composition, dtype={
                "Phase": str,
                "Phase Quantity (wt-%)": np.float64
            })
            compounds = [
                Compound(compound_name,
                    peaks.loc[peaks["Phase"] == compound_name]["Angle (°2?)"],
                    [intensity[(np.abs(angle - peak_angle)).argmin()]
                         for peak_angle in
                         peaks.loc[peaks["Phase"] == compound_name]["Angle (°2?)"]],
                    qtys.loc[qtys["Phase"] == compound_name]["Phase Quantity (wt-%)"].iloc[0])
                for compound_name in compound_names]

        return cls(angle, intensity, compounds)

    @classmethod
    def from_files(cls, ascii_filename, excel_filename, legacy=False):
        if legacy:
            theta2, intensity, _ = read_ascii(ascii_filename)
        else:
            theta2, intensity = read_ascii(ascii_filename)

        meta = pd.read_excel(excel_filename, sheet_name=None)

        compounds = []

        for _, row in meta["info"][["Chemical Formula", "Ref. Code"]].iterrows():
            name = row["Chemical Formula"]
            ref = row["Ref. Code"]

            sample_peaks = meta["sample"]
            matched_peak_positions_raw = np.array(
                sample_peaks[sample_peaks["Matched by"].str.contains(ref, na=False)]["Pos. [°2Th.]"]
            )
            if legacy:
                matched_peak_positions = matched_peak_positions_raw
            else:
                regex = [re.match(r"(\d+\.\d+)\(\d+\)", str(peak))
                    for peak in matched_peak_positions_raw]
                if None in regex:
                    raise Exception("SHIT")
                # it works i promise
                matched_peak_positions = np.array(
                    [float(peak.group(1)) for peak in regex]) # type: ignore

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

    def rename_compound(self, current, new):
        for compound in self.compounds:
            if compound.name == current:
                compound.name = new
                break

        return self


    def subtract_bg(self, filename):

        bg_theta2, bg_intensity = read_ascii(filename)
        # close enough!
        if (
            round(np.min(self.theta2)) != round(np.min(bg_theta2))
            or round(np.max(self.theta2)) != round(np.max(bg_theta2))
            or len(bg_theta2) != len(self.theta2)
        ):
            raise Exception("hmm")

        processed_intensity = self.intensity - bg_intensity
        processed_compounds = [
            compound.subtract_bg(bg_theta2, bg_intensity) for compound in self.compounds]
        

        return XRDSample(self.theta2, processed_intensity, processed_compounds)


    def normalize(self):
        max = np.max(self.intensity)
        normalized_intensity = self.intensity / max
        normalized_compounds = [compound.normalize(max) for compound in self.compounds]

        return XRDSample(self.theta2, normalized_intensity, normalized_compounds)

    def offset(self, amount):
        offset_intensity = self.intensity + amount

        offset_compounds = [compound.offset(amount) for compound in self.compounds]

        return XRDSample(self.theta2, offset_intensity, offset_compounds)

    def plot(self,
             title=None,
             yunits=False,
             exclude_compounds: list = [],
             include_compounds: list = [],
             phase_labels: PhaseLabel = PhaseLabel.LINES,
             n_peaks: int = 8
    ) -> Tuple[Figure, Axes]:
        if ((len(exclude_compounds) + len(include_compounds))
                > abs(len(exclude_compounds) - len(include_compounds))):
            raise Exception("Please only use include or exclude, not both.")

        fig, ax = plt.subplots()
        ax.set_ylabel("Intensity (A.U.)")
        if not yunits:
            ax.set_yticks([])
        ax.set_xlabel("2$\\theta$ (degrees)")
        line = ax.plot(self.theta2, self.intensity, label=title)

        compounds = {
                compound.name: {
                "theta2": compound.peak_theta2, "intensity": compound.peak_intensity}
            for compound in self.compounds
        }

        if len(include_compounds) > 0:
            for key in list(compounds.keys()):
                if key not in include_compounds:
                    del compounds[key]
        match phase_labels:
            case PhaseLabel.MARKERS:
                scatters = xrd_multiplot_label_markers(
                    ax, compounds, abs_peak_offset, margin, n_peaks)

                ax.legend(handles=scatters, loc="upper left", frameon=True)
            case PhaseLabel.LINES:
                lines = xrd_multiplot_label_lines(
                    ax, compounds, n_peaks)
                handles = [line[0] for line in lines]
                ax.legend(handles=handles, loc="upper left", frameon=True)

        return fig, ax

    def max_intensity(self):
        return np.max(self.intensity)


def xrd_multiplot(
    samples: dict[str, XRDSample],
    title: str | None = None,
    margin: float = 0.1,
    peak_label_offset: float = 0.01,
    peak_label_gap: float = 0.07,
    exclude_compounds: list = [],
    include_compounds: list = [],
    phase_labels: PhaseLabel = PhaseLabel.MARKERS,
    n_peaks: int = 8
) -> Tuple[Figure, Axes]:
    if (
        (len(exclude_compounds) + len(include_compounds))
            > abs(len(exclude_compounds) - len(include_compounds))):
        raise Exception("Please only use include or exclude, not both.")

    spacing = 1 + margin

    fig, ax = plt.subplots()#figsize=(12 * 0.6, 8 * 0.6))
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
            color=LINE_COLOURS[i % len(LINE_COLOURS)],
            linewidth=1.3
        )
        for i, (name, sample) in enumerate(offset_samples.items())
    ]

    ax.add_artist(
        ax.legend(
            handles=[item for sublist in lines for item in sublist],
            loc="upper right",
            frameon=True
        )
    )

    compounds = xrd_collate_compounds(offset_samples)
    
    if len(include_compounds) > 0:
        for key in list(compounds.keys()):
            if key not in include_compounds:
                del compounds[key]

    match phase_labels:
        case PhaseLabel.MARKERS:
            scatters = xrd_multiplot_label_markers(
                ax, compounds, abs_peak_offset, margin, n_peaks)

            ax.legend(handles=scatters, loc="center right", frameon=True)
        case PhaseLabel.LINES:
            lines = xrd_multiplot_label_lines(
                ax, compounds, n_peaks)
            handles = [line[0] for line in lines]
            ax.legend(handles=handles, loc="center right", frameon=True)

    return fig, ax


def xrd_collate_compounds(samples: dict[str, XRDSample]) -> dict:
    # collating compounds
    # hate this code
    compounds = {}
    for sample in samples.values():
        for compound in sample.compounds:
            if compound.name not in compounds.keys():
                compounds[compound.name] = {"theta2": [], "intensity": []}
            compounds[compound.name]["theta2"] = np.append(
                compounds[compound.name]["theta2"], compound.peak_theta2
            )
            compounds[compound.name]["intensity"] = np.append(
                compounds[compound.name]["intensity"], compound.peak_intensity
            )
    return compounds


def xrd_multiplot_label_markers(ax, compounds, abs_peak_offset, margin, n_peaks=5):
    # plotting compounds
    # there's probably a pretty functional way of doing this. I do not know it
    scatter = []
    ive_lost_the_plot = []
    for i, (name, value) in enumerate(compounds.items()):
        for angle in value["theta2"][:n_peaks]:
            ive_lost_the_plot.append(angle)

        intensities_offest = [
            intensity + margin * ive_lost_the_plot.count(angle)
            for angle, intensity in zip(value["theta2"][:n_peaks], value["intensity"][:n_peaks])]
        scatter.append(ax.scatter(
            value["theta2"][:n_peaks],
            intensities_offest + abs_peak_offset,
            label=f"\\ce{{{name.replace(" ", "")}}}",
            color="black",
            marker=MARKERS[i],
            zorder=2,
        ))

    return scatter


def xrd_multiplot_label_lines(ax, compounds, n_peaks=5):
    lines = []
    styles = ["solid", "dashed", "dashdot", "dotted"]
    for i, (name, value) in enumerate(compounds.items()):
        lines.append([ax.axvline(
            angle,
            label=f"\\ce{{{name.replace(" ", "")}}}",
            color="grey",
            linestyle=styles[i % 4],
            zorder=-2,
        ) for angle in value["theta2"][:n_peaks]])

    return lines

def xrd_multibar(
    samples: dict[str, XRDSample],
    title: str | None = None,
    exclude_compounds: list = [],
    include_compounds: list = [],
) -> Tuple[Figure, Axes]:
    if (
        (len(exclude_compounds) + len(include_compounds))
            > abs(len(exclude_compounds) - len(include_compounds))):
        raise Exception("Please only use include or exclude, not both.")

    fig, ax = plt.subplots()#figsize=(12 * 0.6, 8 * 0.6))
    if title is not None:
        ax.set_title(title)
    ax.set_ylabel("Phase Composition (Weight \\%)")
    ax.set_xlabel("Sample")

    # tabular data: each column is a sample. each row a compound
    # fill in zeros for samples without data for that compound

    # O(3n) - could probably do all this in one iteration but i think this will
    # be easier to understand.
    # iterate over samples. ideantify all unique compounds
    # iterate over samples again. create qtys for each compound fill zeros -
    # sum on duplicates
    # row stack, transpose
    # pandas it to make it easier
    # | #index | Sample A | Sample B | Sample C |
    # | ph 1   | 50%      |  0%      | 50%      |
    # | ph 2   |  0%      | 50%      | 50%      |
    # iterate over table, plot stacked bar
    unique_phases = set([i for g in
        [[compound.name for compound in sample.compounds]
            for _, sample in samples.items()]
        for i in g])

    data_dict = {phase: [] for phase in unique_phases}

    names = []
    for name, sample in samples.items():
        names.append(name)
        for phase in unique_phases:
            qtys = [compound.qty for compound in sample.compounds if compound.name == phase]
            data_dict[phase].append(sum(qtys))

    bottom = np.zeros(len(names))
    for i, (phase, qtys) in enumerate(data_dict.items()):
        b = ax.bar(names, qtys, 0.6, label=phase, bottom=bottom, color=LINE_COLOURS[i % len(LINE_COLOURS)])
        ax.bar_label(b, label_type="center", fmt=lambda x: f"\\ce{{{phase}}}\n{x:.0f}\\%" if x > 0 else "", color="white", size="x-large")

        bottom += qtys


    return fig, ax

