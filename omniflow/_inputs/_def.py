from __future__ import annotations

import re
import typing
from typing import Callable
import functools

import pypath_common._misc as _common

"""
Define, access and process fields from network inputs.
"""


def _to_fields(loc: dict) -> dict[str, Field]:

    return {k: Field(v) for k, v in loc.items() if k != 'cls'}


class NodeDef(typing.NamedTuple):
    id: Field
    id_type: Field
    entity_type: Field
    organism: Field

    def __new__(
            cls,
            id: Field,
            id_type: Field,
            entity_type: Field,
            organism: Field,
        ):

        args = _to_fields(locals())

        return super().__new__(cls, **args)


class EdgeDef(typing.NamedTuple):
    type: Field
    directed: Field
    effect: Field

    def __new__(cls, type, directed, effect):

        args = _to_fields(locals())

        return super().__new__(cls, **args)


class Field:

    _SEP0 = '?'
    _SEP1 = '!'
    _DEFVALS = (True, False, None)
    _SPECVAL = {'True': True, 'False': False, 'None': None}
    _RECOLPFX = re.compile(r'^__')

    def __init__(
            self,
            definition: str | Callable | tuple | dict,
        ):
        """
        Definition of a field in a data frame.

        Args:
            definition:
                Defines how to extract and process the value of the field from
                a data frame row (record). Most often it refers to a column,
                and defines additional processing steps. The following options
                are supported:
                - A string starting with "__": refers to a column, and it
                  can include further processing by using a simple
                  notation, see more at the ``short`` method of this class.
                - A string: a constant value (the same value will be
                  returned for all any record).
                - A string refering a method in a module, e.g.
                  "json.loads": this method will be used to process the
                  field. By default the function is called with complete
                  records, to operate on a single field instead, add the
                  column name in front, e.g. "__extra_attrs?json.loads".
                - Custom Python code: the code should evaluate into a
                  callable which then will be used to process the field.
                - More or less identical to the above, the definition can
                  be provided as a tuple.
                - A dict representation is also accepted, the dict should
                  contain the column name, the within field separator and
                  the mappings under the "col", "sep" and "maps" keys. If
                  "maps" is callable, it will be called on the field value.
        """

        self._setup(definition)


    def _setup(self, definition: str | Callable | tuple | dict):

        if isinstance(definition, str):

            if definition.startswith('__'):

                # special notation __ marks the string as column definition
                definition = (definition[2:],)

            elif re.match(r'^\w+$', definition):

                # directly provide a constant value
                self._access = lambda record: definition

            elif re.match(r'^[\.\w]+$', definition):

                # refer to any custom function from any module
                self._access = self._func_from_module(definition)

            else:

                # use custom python code
                definition = _common.code_to_func(definition)

        elif isinstance(definition, tuple):

            definition = self.short(self._SEP0.join(definition))

        if isinstance(definition, Callable):

            self._access = definition

        elif isinstance(definition, dict):

            self._access = functools.partial(
                self._access_by_param,
                param = definition,
            )


    @classmethod
    def short(cls, definition: str) -> dict:
        """
        Parse field processing parameters from a short string notation.

        The simple notation is specific to this module and makes it possible to
        easily define the most common, basic ways of processing data from
        columns. For the simplest cases not even this notation is required, for
        more complex cases you can always inject your custome Python code.

        The notation uses `?` and `!` as higher and a lower level separators,
        hence these characters can not be part of column names or separators
        used in the processed data or in the definition. The first part of the
        notation describes the column name, the separators used to split the
        column values and to split the values in the definition. The rest of
        the definition consists of potential mappings. The exact processing
        behaviour depends on the number and contents of the mappings:

        __column = Simply take the value of "column".
        __column!; = Split the value of "column" by ";".
        __column!;?Activation = Split the value of "column" and if the word
            "Activation" is among the values return True, otherwise False.
        __column!!|?Activation?Inhibition?Undefined|Unknown = Do not split the
            value of column, but split the values in this definition by "|".
            Then match the column value against the ones listed here, i.e.
            Activation = True, Inhibition = False, Undefined or Unknown = None.
            These (True, False, None) in this order are values of the default
            mapping, different values have to be specified explicitely.
        __column!!,?Ser!Serine,S,Ser?Tyr!Tyrosine,Tyr,Y?Thr,Threonine,Thr,T =
            Map single and three letter notation and full names of these amino
            acides to their three letter notation.

        In this function the "__" prefix is optional, as here we can be sure we
        refer to columns.
        In rare cases when the data is incompatible with `?` and `!` being used
        as separators, change the ``_SEP0`` and the ``_SEP1`` attributes of
        this class.
        """

        definition = _RECOLPFX.sub('', definition)
        col, *maps = definition.split(cls._SEP0)
        col, *seps = col.split(cls._SEP1)
        datasep, defsep = (s or None for s in (seps + [None] * 2)[:2])

        maps =
            self._func_from_module(maps)
                if maps and re.match(r'^[\w\.]+$', maps[0])
            [self._process_map(i, m, defsep) for i, m in enumerate(maps)]
        )

        return {
            'col': col,
            'sep': datasep,
            'maps': maps,
        }


    @classmethod
    def _process_map(
            cls,
            i: int,
            m: str,
            defsep: str | None = None,
        ) -> tuple[str | set[str], Any]:

        values, *map_to = m.split(cls._SEP1)[::-1] + [None]
        map_to = map_to[0] or cls._DEFVALS[i]
        map_to = cls._SPECVAL.get(map_to, map_to)
        values = set(values.split(defsep)) if defsep else values

        return values, map_to

    @classmethod
    def _access_by_param(
            cls,
            record: pd.DataFrame | dict,
            param: str | dict,
        ) -> Any:

        param = cls.short(param) if isinstance(param, str) else param
        return cls.access_by_param(record, param)


    @staticmethod
    def access_by_param(record: pd.DataFrame | dict, param: dict) -> Any:
        """
        Access the field value using processing parameters.

        Args:
            record:
                One row from a data frame or a dict.
            param:
                Processing parameters: a dict with keys "col", "sep" and
                "maps".

        Returns:
            The processed value of the field.
        """

        value = record[param['col']]
        op_one, op_many = value.__eq__, lambda x: x.__contains__(value)

        if param['sep']:

            value = set(value.split(param['sep']))
            op_one, op_many = value.__contains__, value.__and__

        if callable(param['maps']):

            return param['maps'](value)

        for v, map_to in param['maps']:

            if match := (op_many(v) if isinstance(v, set) else op_one(v)):

                return map_to

        return value


    @staticmethod
    def _access_column(col: str, record: pd.DataFrame | dict) -> Any:

        return record[col]


    @staticmethod
    def _func_from_module(func: str) -> Callable:

        # refer to any custom function from any module
        mod, *func = definition.rsplit('.', maxsplit = 1)
        mod = __import__(mod)
        return getattr(mod, func)


    def __call__(self, record: pd.DataFrame | dict) -> Any:

        return self._access(record)
