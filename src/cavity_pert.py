import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd
from dataclasses import dataclass
from typing import Self

from pandas.io.formats.format import math


def real_perm(g, vc, vs, f0, fs):
    return (2 * ((f0 - fs) / f0) * (vc / vs) * g) + 1


def imag_perm(g, vc, vs, q0, qs):
    return ((1 / qs) - (1 / q0)) * (vc / vs) * g


def temperature_correct_f(temp, f_sample, f_correction_mode):
    # vastly different results depending on temperature used for correct: I'm
    # doing it based on the minimum, which should usually be about 25 degrees.
    min_index = np.argmin(temp)
    return f_sample / (f_correction_mode / f_correction_mode[min_index])

@dataclass
class CavityData:
    radius: float 
    height: float
    modes: dict[str, int]

    def volume(self):
        return self.height*math.pi*self.radius**2

    def sample_volume(self, sample_radius):
        return self.height*math.pi*(sample_radius**2)

@dataclass
class Cavities:
    # eps-cav, from Jerome's thesis.
    GANDALF = CavityData(
        radius=0.046,
        height=0.04,
        modes={"TM010": 0.294,
        "TM011": 0.962,
        "TM020": 0.124,
        "TM021": 0.183,
        "TM012": 2.580,
        "TM022": 0.332}
    )

    # mu-cav, from Jerome's thesis.
    SARUMAN = CavityData(
        radius=0.0475,
        height=0.04,
        modes={"TM110": 0.388,
        "TE011": 0.285,
        "TE021": 0.107,
        "TE012": 0.494,
        "TE022": 0.130}
    )

    Clamshell2495Mag = CavityData(
        radius=0.0460,
        height=0.04,
        modes={"TM010": 0.889365948}
    )


class PermResult:
    def __init__(self, data, fn, coeff) -> None:
        self.data = data
        self.fn = fn
        self.coeff = coeff
        self.mean = np.mean(data)
        self.max = np.max(data)
        self.min = np.min(data)

@dataclass
class ModeData:
    name: str
    frequency: np.typing.ArrayLike
    q: np.typing.ArrayLike

@dataclass
class CavityPerturbationMeasurement:
    measurement_modes: list[ModeData]
    temperature_modes: list[ModeData]
    temperature: np.ndarray
    frequency: np.ndarray

    @classmethod
    def from_tdms(cls, filepath) -> Self:
        return cls()

    @classmethod
    def from_csv(cls, filepath) -> Self:
        df = pd.read_csv(filepath)
        # extract modes from column headers
        f_regex = re.compile(r"Frequency ([\w\d\-\_]*)\/Hz")
        regexed = [f_regex.match(col) for col in df.columns]
        modes = [result.group(1) for result in regexed if result]
        print(modes)

        # bisect on suffix to find temperature modes
        standard_mode_names = [mode for mode in modes if mode[-1:-2] != "_t"]
        temp_mode_names = [mode for mode in modes if mode[-1:-2] == "_t"]
        print(standard_mode_names)
        print(temp_mode_names)

        # extract frequency series
        # extract q series

        return cls()


class CavityPerturbationResult:
    def __init__(self, real: PermResult, imag: PermResult, temperature, mean_f,
                 f_shift, q_shift) -> None:
        self.imag = imag
        self.real = real
        self.temperature = temperature
        self.mean_f = mean_f
        self.f_shift = f_shift
        self.q_shift = q_shift

    def plot(self, title):
        _, ax = plt.subplots()
        ax.set_title(title)
        ax.set_xlabel("Temperature, Celcius")

        ax.plot(self.temperature, self.real.data, label="Real component, raw.")
        ax.plot(
            self.temperature,
            self.real.fn(self.temperature),
            label="Real component, fit.",
        )
        ax.set_ylabel("Real component")
        ax.legend()

        real_min = self.real.data.min()
        real_max = self.real.data.max()
        real_range = real_max - real_min

        ax.set_ylim(real_min - real_range, real_max)

        imag_ax = ax.twinx()
        imag_ax.plot(self.temperature, self.imag.data, color="red", label="Imag component, raw.")  # type: ignore
        imag_ax.plot(self.temperature, self.imag.fn(self.temperature), color="green", label="Imag component, fit.")  # type: ignore
        imag_ax.legend()
        imag_ax.set_ylabel("Imaginary component")

        imag_min = self.imag.data.min()
        imag_max = self.imag.data.max()
        imag_range = imag_max - imag_min

        imag_ax.set_ylim(imag_min, imag_max + imag_range)


def analyse(
    sample_tdms_file, empty_tdms_file, cavity, sample_radius, temp_correct=True, use_modes=None
): #-> dict[str, CavityPerturbationResult]:

    sample_modes = []
    temp_correct_modes = []
    temperature = None
    for group in sample_tdms_file.groups():
        if group.name == "Temperature":
            temperature = group["Ch0"][:]
            continue

        if re.search(".*_t", group.name):
            temp_correct_modes.append(group)
        else:
            if use_modes is not None:
                if group.name in use_modes:
                    sample_modes.append(group)
            else:
                sample_modes.append(group)

    if temperature is None:
        raise Exception("No temp series in TDMS file")


    results = {}
    input_data = {}
    for mode in sample_modes:
        empty_f0 = np.average(empty_tdms_file[mode.name]["F0 (Hz)"][:])
        empty_q0 = np.average(empty_tdms_file[mode.name]["Q0"][:])

        f = mode["F0 (Hz)"][:]
        q = mode["Q0"][:]
        mode_volume = cavity.modes[mode.name]
        cavity_volume = cavity.volume()
        sample_volume = cavity.sample_volume(sample_radius)
        print(f"{mode.name=}, {mode_volume=}, {cavity_volume=}, {sample_volume=}")

        if temp_correct:
            f_correction_mode = temp_correct_modes[0]["F0 (Hz)"][:]
            f_corrected = temperature_correct_f(temperature, f, f_correction_mode)

        sample_real_perm = real_perm(
            mode_volume, cavity_volume, sample_volume, empty_f0,
            f_corrected if temp_correct else f)

        sample_imag_perm = imag_perm(
            mode_volume, cavity_volume, sample_volume, empty_q0, q
        )

        real_perm_coeff = np.polyfit(temperature, sample_real_perm, 1)
        imag_perm_coeff = np.polyfit(temperature, sample_imag_perm, 1)
        real_perm_fn = np.poly1d(real_perm_coeff)
        imag_perm_fn = np.poly1d(imag_perm_coeff)

        results[mode.name] = CavityPerturbationResult(
            PermResult(sample_real_perm, real_perm_fn, real_perm_coeff),
            PermResult(sample_imag_perm, imag_perm_fn, imag_perm_coeff),
            temperature,
            np.mean(f),
            (empty_f0 - f) / empty_f0,
            (1 / q) - (1 / empty_q0),
        )
        input_data[mode.name] = {
            "f": f,
            "q": q
        }

    return (results, input_data)
