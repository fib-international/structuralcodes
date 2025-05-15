"""A concrete implementation of a shell geometry."""

import typing as t

import numpy as np

from ..core.base import Material
from ._geometry import Geometry


class ShellReinforcement(Geometry):
    """A class for representing reinforcement in a shell geometry."""

    _z: float
    _n_bars: float
    _cc_bars: float
    _diameter_bar: float
    _material: Material
    _phi: float

    def __init__(
        self,
        z: float,
        n_bars: float,
        cc_bars: float,
        diameter_bar: float,
        material: Material,
        phi: float,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
    ) -> None:
        """Initialize a shell reinforcement."""
        super().__init__(name, group_label)

        self._z = z
        self._n_bars = n_bars
        self._cc_bars = cc_bars
        self._diameter_bar = diameter_bar
        self._material = material
        self._phi = phi

    @property
    def z(self) -> float:
        """Return the reinforcement position over the thickness."""
        return self._z

    @property
    def n_bars(self) -> float:
        """Return the number of bars per unit width."""
        return self._n_bars

    @property
    def cc_bars(self) -> float:
        """Return the spacing between bars."""
        return self._cc_bars

    @property
    def diameter_bar(self) -> float:
        """Return the bar diameter."""
        return self._diameter_bar

    @property
    def material(self) -> Material:
        """Return the material of the reinforcement."""
        return self._material

    @property
    def phi(self) -> float:
        """Return the orientation angle of the reinforcement."""
        return self._phi

    @property
    def T(self) -> np.ndarray:
        """Return the transformation matrix for the reinforcement."""
        c, s = np.cos(self.phi), np.sin(self.phi)
        return np.array(
            [
                [c * c, s * s, c * s],
                [s * s, c * c, -c * s],
                [-2 * c * s, 2 * c * s, c * c - s * s],
            ]
        )

    def _repr_svg_(self) -> str:
        """Returns the svg representation."""
        raise NotImplementedError


class ShellGeometry(Geometry):
    """A class for a shell with a thickness and material."""

    _reinforcement: t.List[ShellReinforcement]

    def __init__(
        self,
        thickness: float,
        material: Material,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
    ) -> None:
        """Initialize a shell geometry."""
        super().__init__(name=name, group_label=group_label)

        if thickness <= 0:
            raise ValueError('Shell thickness must be positive.')

        self._thickness = thickness
        self._material = material

        self._reinforcement = []

    @property
    def thickness(self) -> float:
        """Return the shell thickness."""
        return self._thickness

    @property
    def material(self) -> Material:
        """Return the material of the shell."""
        return self._material

    @property
    def reinforcement(self) -> t.List[ShellReinforcement]:
        """Return all reinforcement layers."""
        return self._reinforcement

    def add_reinforcement(
        self,
        reinforcement: t.Union[ShellReinforcement, t.List[ShellReinforcement]],
    ) -> None:
        """Add reinforcement to the shell geometry."""
        if isinstance(reinforcement, ShellReinforcement):
            self._reinforcement.append(reinforcement)
            to_validate = [reinforcement]
        elif isinstance(reinforcement, list):
            self._reinforcement.extend(reinforcement)
            to_validate = reinforcement
        else:
            raise TypeError(
                'Reinforcement must be a ShellReinforcement or a list of '
                'ShellReinforcement.'
            )

        # Validate each reinforcement layer
        for r in to_validate:
            half_thickness = self._thickness / 2
            if not (
                -half_thickness + r.diameter_bar / 2
                <= r.z
                <= half_thickness - r.diameter_bar / 2
            ):
                raise ValueError(
                    f'Reinforcement at z = {r.z:.2f} mm is outside the '
                    f'range [-{half_thickness:.2f}, {half_thickness:.2f}] mm.'
                )

    def _repr_svg_(self) -> str:  # noqa: PLR0915
        """Returns the svg representation."""
        # overall drawing area dimensions
        w, h = 1000, 1000
        t_val = self.thickness
        # half-width/height for centering
        hw, hh = w / 2, h / 2
        # single view dimensions
        sw, sh = w + 200, h + 200
        # side view height based on thickness
        sh_side = int(t_val + 200)
        # main and side view ViewBox definitions
        vb_main = f'{-hw - 100} {-hh - 100} {w + 200} {h + 200}'
        vb_side = f'{-hw - 100} {-t_val / 2 - 100} {w + 200} {t_val + 200}'
        # scaling and spacing between views
        scale, gap = 0.6, 50

        def draw_rebars(view: str) -> str:
            # draw parallel rebar lines for top/bottom plan views
            lines: t.List[str] = []
            for layer in self.reinforcement:
                # select only top or bottom layers
                top = view == 'top' and layer.z >= 0
                bot = view == 'bottom' and layer.z < 0
                if not (top or bot):
                    continue
                # orientation and spacing
                phi, sp = layer.phi, layer.cc_bars
                c, s = np.cos(phi), np.sin(phi)
                # perpendicular offset vector
                px, py = -s, c
                # color by orientation
                col = (
                    'red'
                    if np.isclose(phi % np.pi, 0)
                    else 'blue'
                    if np.isclose(phi % np.pi, np.pi / 2)
                    else 'green'
                )
                n = int(w / sp) + 3  # linecount with safety margin to angles
                L = 2000  # line length extension
                for i in range(-n // 2, n // 2 + 1):
                    ox, oy = i * sp * px, i * sp * py
                    # line endpoints
                    x1, y1 = ox - L * c, oy - L * s
                    x2, y2 = ox + L * c, oy + L * s
                    for _ in range(int(layer.n_bars)):
                        # actual bar offset
                        lines.append(
                            f'<line x1="{x1}" y1="{y1}" '
                            f'x2="{x2}" y2="{y2}" '
                            f'stroke="{col}" '
                            f'stroke-width="{layer.diameter_bar}" '
                            'stroke-opacity="0.7"/>'
                        )
            return ''.join(lines)

        def build_plan(view: str) -> str:
            # define clipping area for concrete
            clip = (
                '<defs><clipPath id="clip' + view + '">'
                f'<rect x="{-hw}" y="{-hh}" '
                f'width="{w}" height="{h}"/>'
                '</clipPath></defs>'
            )
            # concrete background
            bg = (
                f'<rect x="{-hw}" y="{-hh}" '
                f'width="{w}" height="{h}" '
                'fill="lightgray"/>'
            )
            # concrete outline
            out = (
                f'<rect x="{-hw}" y="{-hh}" '
                f'width="{w}" height="{h}" '
                'fill="none" stroke="black" '
                'stroke-width="2"/>'
            )
            # view label
            lbl = (
                f'<text x="{-hw + 20}" y="{-hh - 40}" '
                'font-size="30">' + f'{view} view</text>'
            )
            # rebar lines clipped to concrete area
            body = (
                '<g clip-path="url(#clip'
                + view
                + ')">'
                + draw_rebars(view)
                + '</g>'
            )
            return clip + bg + body + out + lbl

        def build_side(ax: str) -> str:
            # side view perpendicular direction
            span = w if ax == 'x' else h
            half = span / 2
            t_ht = self.thickness
            cid = 'clipSide' + ax
            # clipping region for shell thickness
            clip = (
                '<defs><clipPath id="' + cid + '">'
                f'<rect x="{-half}" y="{-t_ht / 2}" '
                f'width="{span}" height="{t_ht}"/>'
                '</clipPath></defs>'
            )
            # shell background
            bg = (
                f'<rect x="{-half}" y="{-t_ht / 2}" '
                f'width="{span}" height="{t_ht}" '
                'fill="lightgray"/>'
            )
            # shell outline
            out = (
                f'<rect x="{-half}" y="{-t_ht / 2}" '
                f'width="{span}" height="{t_ht}" '
                'fill="none" stroke="black" '
                'stroke-width="2"/>'
            )
            # side view label
            lbl = (
                f'<text x="{-half + 20}" y="{-t_ht / 2 - 40}" '
                'font-size="30">' + f'Side {ax.upper()}-view</text>'
            )
            elems: t.List[str] = []
            for layer in self.reinforcement:
                # select bars not parallel to this side-plane
                perp = not (
                    (
                        ax == 'x'
                        and np.isclose(layer.phi % np.pi, np.pi / 2, atol=1e-2)
                    )
                    or (
                        ax == 'y'
                        and np.isclose(layer.phi % np.pi, 0, atol=1e-2)
                    )
                )
                if perp:
                    # draw circles for perp layers
                    step = layer.cc_bars or span
                    kmax = int(np.floor(half / step))
                    xs = [k * step for k in range(-kmax, kmax + 1)]
                    col = (
                        'red'
                        if np.isclose(layer.phi % np.pi, 0)
                        else 'blue'
                        if np.isclose(layer.phi % np.pi, np.pi / 2)
                        else 'green'
                    )
                    for x in xs:
                        for j in range(int(layer.n_bars)):
                            ex = (
                                j - (layer.n_bars - 1) / 2
                            ) * layer.diameter_bar
                            cy = -layer.z
                            elems.append(
                                f'<circle cx="{x + ex:.1f}" '  # bar position
                                f'cy="{cy:.1f}" '
                                f'r="{layer.diameter_bar / 2}" '
                                f'fill="{col}" '
                                'fill-opacity="0.9"/>'
                            )
                else:
                    # draw lines for parallel layers
                    col = (
                        'red'
                        if np.isclose(layer.phi % np.pi, 0)
                        else 'blue'
                        if np.isclose(layer.phi % np.pi, np.pi / 2)
                        else 'green'
                    )
                    y = -layer.z
                    for j in range(int(layer.n_bars)):
                        ex = (j - (layer.n_bars - 1) / 2) * layer.diameter_bar
                        x1 = -half + ex
                        x2 = half + ex
                        elems.append(
                            f'<line x1="{x1}" y1="{y}" '  # bar top view
                            f'x2="{x2}" y2="{y}" '
                            f'stroke="{col}" '
                            f'stroke-width="{layer.diameter_bar}" '
                            'stroke-opacity="0.7"/>'
                        )
            body = (
                '<g clip-path="url(#' + cid + ')">' + ''.join(elems) + '</g>'
            )
            return clip + bg + body + out + lbl

        # assemble all SVG fragments
        svg = [
            f'<svg width="{(sw * 2 + gap) * scale:.0f}" '
            f'height="{(sh + gap + sh_side) * scale:.0f}" '
            f'viewBox="0 0 {(sw * 2 + gap):.0f} '
            f'{(sh + gap + sh_side):.0f}"'
            ' xmlns="http://www.w3.org/2000/svg"><g>'
        ]
        # plan views
        svg.append(
            '<g><svg width="' + str(sw) + '" '  # top plan
            'height="' + str(sh) + '" '
            'viewBox="' + vb_main + '">' + build_plan('top') + '</svg></g>'
        )
        svg.append(
            '<g transform="translate(' + str(sw + gap) + ',0)">'  # bottom plan
            '<svg width="' + str(sw) + '" '
            'height="' + str(sh) + '" '
            'viewBox="' + vb_main + '">' + build_plan('bottom') + '</svg></g>'
        )
        # side views
        svg.append(
            '<g transform="translate(0,' + str(sh + gap) + ')">'  # side-x
            '<svg width="' + str(sw) + '" '
            'height="' + str(sh_side) + '" '
            'viewBox="' + vb_side + '">' + build_side('x') + '</svg></g>'
        )
        svg.append(
            '<g transform="translate('
            + str(sw + gap)
            + ','
            + str(sh + gap)
            + ')">'  # side-y
            '<svg width="' + str(sw) + '" '
            'height="' + str(sh_side) + '" '
            'viewBox="' + vb_side + '">' + build_side('y') + '</svg></g>'
        )
        svg.append('</g></svg>')
        return ''.join(svg)
