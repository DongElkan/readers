import re
import base64
import zlib
import collections
import numpy as np

from operator import attrgetter

from .base import MassSpectrum, XIC

try:
    import xml.etree.cElementTree as et
except ImportError:
    import xml.etree.ElementTree as et


def _decode_binary(peak_string, precision, is_decompressed):
    """
    Decode binary string to float points.
    """
    # set up data type
    float_ = np.float32 if precision == 32 else np.float64
    dtype = np.dtype([('m/z array', float_),
                      ('intensity array', float_)]).newbyteorder('>')

    # decode peak string
    decoded = base64.b64decode(peak_string)
    if is_decompressed:
        decoded = zlib.decompress(decoded)

    peak_conv = np.frombuffer(bytearray(decoded), dtype=dtype)

    # numpy array peak list
    return np.array(list(zip(*peak_conv.tolist()))).T


class mzXMLReader():
    """ Parse mzXML files to read mass spectral peaks. """
    def __init__(self, mzxml_file):
        self.mzxml_file = mzxml_file
        # read mass spectra
        self._read_spectrum()
        # organize them based on ms_levels
        self._ms_levels = collections.defaultdict()
        for level in range(1, 5):
            idx = set(spectrum.scan for spectrum in self._spectra
                      if spectrum.ms_level == level)
            if not idx:
                break
            self._ms_levels[level] = idx

    def get_mass_spectrum(self, index):
        """
        Get tandem mass spectra according to the input scans.
        If the no scan is specified, all tandem mass spectra
        are output.

        Parameters
        ----------
        index: None, array or list
               Index of spectra extracted. If is set to None,
               all mass spectra are output.
        """
        if index is None:
            return self._spectra

        if isinstance(index, int):
            return self._spectra[index]

        return [self._spectra[i] for i in sorted(index)]

    def get_xic(self, mz, tol=0.1, unit="Da"):
        """
        Get MS1 according to the input scans

        Parameters
        ----------
        mz: array or list
            m/z values for extracting XIC.
        tol: float
             Tolerance for extracting XIC.
        unit: Da | ppm
              Unit of tolerance, default is "Da"

        Returns
        -------
        XIC objects

        """
        t = mz / 1e6 * tol if unit == "ppm" else tol
        return self._xic(mz, t)

    def _read_spectrum(self):
        """
        Get spectrum according to input scans
        """
        spectra = []
        for _, elem in et.iterparse(self.mzxml_file):
            if elem.tag.endswith("peaks"):
                peaks = self._parse_spectrum_peaks(elem)

            elif elem.tag.endswith("precursorMz"):
                precursor_mz = float(elem.text)
                charge = elem.get("precursorCharge")
                if charge is not None:
                    charge = int(charge)
                collision = elem.get("activationMethod")

            elif elem.tag.endswith("scan"):
                spectrum = self._get_spectrum_info(elem)
                spectrum.spectrum = peaks
                # assign precursors and charge to MS2
                if spectrum.ms_level == 2:
                    spectrum.charge = charge
                    spectrum.precursormz = precursor_mz
                    if collision is not None:
                        spectrum.collision = collision
                spectra.append(spectrum)

        self._spectra = sorted(spectra, key=attrgetter("scan"))

    def _get_spectrum_info(self, scan_info):
        """ Get spectral information. """
        rt = float(re.findall(r"PT(.+)S", scan_info.get("retentionTime"))[0])
        mslevel = int(scan_info.get("msLevel"))
        energy = scan_info.get("collisionEnergy")
        if energy is not None:
            energy = float(energy)
        collision = scan_info.get("activationMethod")
        scan = int(scan_info.get("num"))

        return MassSpectrum(**{
            "scan": scan, "rt": rt, "collision": collision, "energy": energy,
            "ms_level": mslevel
        })

    def _parse_spectrum_peaks(self, xml_peak):
        """ Parse spectral peaks from the XML peak string. """
        if xml_peak.text is None:
            return None

        precision = int(xml_peak.get("precision"))
        is_decompressed = xml_peak.get("compressionType") == "zlib"
        peaks = _decode_binary(xml_peak.text, precision, is_decompressed)
        return peaks[peaks[:, 1] > 0, :]

    def _xic(self, mz, tol):
        """ Generate XIC according to matched peaks """
        xic = []
        for spectrum in self._spectra:
            if spectrum.ms_level != 1:
                continue
            ix = np.absolute(spectrum.spectrum[:, 0] - mz) <= tol
            xic.append([spectrum.rt, spectrum.spectrum[ix, 1].max()]
                       if ix.any() else [spectrum.rt, 0])

        return XIC(mz=mz, tol=tol, xic=np.array(xic))
