import numpy as np
import matplotlib.pyplot as plt
import re


def real_perm(g, vc, vs, f0, fs):
    return (2 * ((f0 - fs) / f0) * (vc / vs) * g) + 1


def imag_perm(g, vc, vs, q0, qs):
    return ((1 / qs) - (1 / q0)) * (vc / vs) * g


def temperature_correct_f(temp, f_sample, f_correction_mode):
    # vastly different results depending on temperature used for correct: I'm
    # doing it based on the minimum, which should usually be about 25 degrees.
    min_index = np.argmin(temp)
    return f_sample / (f_correction_mode / f_correction_mode[min_index])


# Cuenca
MODE_VOLUMES = {
    "TM110": 0.388,
    "TE011": 0.271,
    "TE021": 0.108,
    "TE012": 0.598,
    "TE022": 0.162,
}

CAVITY_SQUARED_RADII = {"Gandalf": 0.046**2}


class PermResult:
    def __init__(self, data, fn, coeff) -> None:
        self.data = data
        self.fn = fn
        self.coeff = coeff


class CavityPerturbationResult:
    def __init__(self, real: PermResult, imag: PermResult, temperature, mean_f) -> None:
        self.imag = imag
        self.real = real
        self.temperature = temperature
        self.mean_f = mean_f

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
    sample_tdms_file, empty_tdms_file, cavity_volume, sample_volume
) -> dict[str, CavityPerturbationResult]:

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
            sample_modes.append(group)

    if temperature is None:
        raise Exception("No temp series in TDMS file")

    results = {}
    for mode in sample_modes:
        empty_f0 = np.average(empty_tdms_file[mode.name]["F0 (Hz)"][:])
        empty_q0 = np.average(empty_tdms_file[mode.name]["Q0"][:])

        f = mode["F0 (Hz)"][:]
        f_correction_mode = temp_correct_modes[0]["F0 (Hz)"][:]
        q = mode["Q0"][:]
        mode_volume = MODE_VOLUMES[mode.name]

        f_corrected = temperature_correct_f(temperature, f, f_correction_mode)

        sample_real_perm = real_perm(
            mode_volume, cavity_volume, sample_volume, empty_f0, f_corrected
        )
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
        )

    return results
