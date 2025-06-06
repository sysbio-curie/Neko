from __future__ import annotations

from typing import Literal


EFFECT_TYPES = {
    "1": "stimulation",
    "activate": "stimulation",
    "stimulate": "stimulation",
    "phosphorylate": "stimulation",
    "stimulation": "stimulation",
    "->": "stimulation",
    "-|": "inhibition",
    "-1": "inhibition",
    "inhibit": "inhibition",
    "block": "inhibition",
    "inhibition": "inhibition",
    "form complex": "form_complex",
    "form_complex": "form_complex",
    "form-complex": "form_complex",
    "complex formation": "form_complex"
}


def effect(
        value: str | int
    ) -> Literal['stimulation', 'inhibition', 'form_complex', 'undefined']:
    """
    Determine the effect based on the interaction type.

    Args:
        interaction_type:
            A string representing the type of interaction.

    Return:
       A string representing the effect of the interaction.
       If the interaction type is not recognized, it returns "undefined".
    """

    return EFFECT_TYPES.get(str(value).lower(), "undefined")
