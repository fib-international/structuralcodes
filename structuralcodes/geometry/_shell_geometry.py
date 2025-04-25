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
        elif isinstance(reinforcement, list):
            self._reinforcement.extend(reinforcement)

        # Validate each reinforcement layer
        for r in reinforcement:
            half_thickness = self._thickness / 2
            if not (
                -half_thickness + r.diameter_bar / 2
                <= r.z
                <= half_thickness - r.diameter_bar / 2
            ):
                raise ValueError(
                    f'Reinforcement at z = {r.z:.2f} mm is outside the'
                    f'range [-{half_thickness:.2f}, {half_thickness:.2f}] mm.'
                )

    def _repr_svg_(self) -> str:
        """Returns the svg representation."""
        # Concrete dimensions and half sizes
        w, h = 1000, 1000
        hw, hh = w / 2, h / 2
        # ViewBox for inner SVG elements
        vb = f'{-hw - 100} {-hh - 100} {w + 200} {h + 200}'

        # Draw extended reinforcement lines
        def draw_rebars(view: str) -> str:
            elems = []
            for layer in self.reinforcement:
                # Choose layers based on z-value (top vs bottom)
                if (view == 'top' and layer.z >= 0) or (
                    view == 'bottom' and layer.z < 0
                ):
                    phi = layer.phi
                    sp = layer.cc_bars
                    c, s = np.cos(phi), np.sin(phi)
                    # Perpendicular vector for shifting parallel lines
                    px, py = -s, c
                    # Color based on orientation
                    col = (
                        'red'
                        if np.isclose(phi % np.pi, 0)
                        else 'blue'
                        if np.isclose(phi % np.pi, np.pi / 2)
                        else 'green'
                    )
                    n = int(w / sp) + 3  # Safety margin to take care of angels
                    L = 2000  # Extend lines
                    for i in range(-n // 2, n // 2 + 1):
                        ox = i * sp * px
                        oy = i * sp * py
                        x0, y0 = ox, oy
                        # Extended endpoints along the rebar direction
                        x1, y1 = x0 - L * c, y0 - L * s
                        x2, y2 = x0 + L * c, y0 + L * s
                        for j in range(int(layer.n_bars)):
                            extra = (
                                j - (layer.n_bars - 1) / 2
                            ) * layer.diameter_bar
                            sx, sy = extra * px, extra * py
                            p1 = (x1 + sx, y1 + sy)
                            p2 = (x2 + sx, y2 + sy)
                            elems.append(
                                f'<line x1="{p1[0]}" y1="{p1[1]}" x2="{p2[0]}"'
                                f'y2="{p2[1]}" stroke="{col}" stroke-width="'
                                f'{layer.diameter_bar}" stroke-opacity="0.7"/>'
                            )
            return ''.join(elems)

        # Build one view (top or bottom)
        def build_view(view: str) -> str:
            # Define a clipPath for the concrete area
            clip_def = (
                f'<defs><clipPath id="clipConcrete">'
                f'<rect x="{-hw}" y="{-hh}" width="{w}" height="{h}" />'
                f'</clipPath></defs>'
            )
            # Draw rebar lines and clip them to the concrete area
            rebar_svg = (
                f'<g clip-path="url(#clipConcrete)">{draw_rebars(view)}</g>'
            )
            # Draw a background rectangle for the concrete (lightgray fill)
            bg_rect = (
                f'<rect x="{-hw}" y="{-hh}" width="{w}" height="{h}" '
                + 'fill="lightgray" />'
            )

            # Draw a concrete outline and a label
            outline = (
                f'<rect x="{-hw}" y="{-hh}" width="{w}" height="{h}" '
                + 'fill="none" stroke="black" stroke-width="2" />'
            )
            lbl = (
                f'<text x="{-hw + 20}" y="{-hh - 40}" font-size="30" '
                + 'fill="black">'
                + f'{"Top View" if view == "top" else "Bottom View"}'
                + '</text>'
            )

            return clip_def + bg_rect + rebar_svg + outline + lbl

        scale = 0.6

        # Assemble both views side by side
        gap = 50  # gap between views
        sw, sh = w + 200, h + 200  # single view width/height
        total_w, total_h = (sw * 2 + gap) * scale, sh * scale

        svg_parts = [
            f'<svg width="{total_w * scale}" height="{total_h * scale}" '
            f'viewBox="0 0 {total_w} {total_h}" '
            'xmlns="http://www.w3.org/2000/svg">'
            # everything is grouped to scale
            f'<g transform="scale({scale})">'
        ]

        # Top-view (left)
        svg_parts.append(
            "<g transform='translate(0,0)'><svg x='0' y='0' "
            f"width='{sw}' height='{sh}' viewBox='{vb}'>"
            f"{build_view('top')}</svg></g>"
        )

        # Bottom-view (right)
        svg_parts.append(
            f"<g transform='translate({sw + gap},0)'>"
            f"<svg x='0' y='0' width='{sw}' height='{sh}' "
            f"viewBox='{vb}'>{build_view('bottom')}</svg></g>"
        )

        svg_parts.append('</g></svg>')
        return ''.join(svg_parts)
