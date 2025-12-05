from __future__ import annotations

import logging
import sys
import typing
from dataclasses import dataclass, field

from .context import get_params
from .proxy import (
    struct_to_proxy_class_name,
)
from .structs import (
    ParsedStructure,
    StructureMember,
)
from .transforms import Transform, get_type_transform
from .types import (
    ALLOC,
    INT8,
    NOT,
    PTR,
    REAL16,
    ArgumentType,
    FullType,
    PointerType,
)

if typing.TYPE_CHECKING:
    from .routines import FortranRoutine


logger = logging.getLogger(__name__)


@dataclass(kw_only=True)
class Argument:
    """
    Represents a Fortran structure member, originally named an Argument in the
    cpp_bmad_interface.

    Attributes
    ----------
    is_component : bool
        Whether this is a structure component. If False, it's an array bound.
    f_name : str
        Fortran side name of argument (lowercase).
    c_name : str
        C++ side name of argument, potentially mangled to avoid reserved word conflicts.
    type : str
        Fortran type without parameters, e.g., 'real', 'type', 'character'.
    kind : str
        Fortran kind, e.g., '', 'rp', 'coord_struct'.
    pointer_type : str
        Pointer type: NOT, PTR, or ALLOC.
    array : List[str]
        Array dimension specifications, e.g., [':', ':'] or ['0:6', '3'].
    lbound : List[Any]
        Lower bounds for each array dimension.
    ubound : List[Any]
        Upper bounds for each array dimension.
    init_value : str
        Initialization value.
    comment : str
        Comment from the Fortran structure definition.
    transform : FortranSideTransform
        Fortran/C++ side type translation.
    """

    is_component: bool = True
    f_name: str = ""
    c_name: str = ""
    type: ArgumentType = "real"
    kind: str = ""
    pointer_type: PointerType = NOT
    array: list[str] = field(default_factory=list)
    init_value: str | None = None
    comment: str = ""
    member: StructureMember

    @property
    def is_pointer(self) -> bool:
        return self.pointer_type == {"PTR", "ALLOC"}

    @property
    def is_dynamic_array(self) -> bool:
        return ":" in self.array or "0:" in self.array or "*" in self.array

    @classmethod
    def from_fstruct(cls, fstruct: ParsedStructure | FortranRoutine, member: StructureMember):
        if member.kind and member.type.lower() == "integer":
            type_ = INT8
        else:
            type_ = member.type

        if member.kind and member.kind.lower() == "qp" and member.type.lower() == "real":
            type_ = REAL16

        if member.type_info.pointer:
            pointer_type = PTR
        elif member.type_info.allocatable:
            pointer_type = ALLOC
        else:
            pointer_type = NOT

        f_to_c_name = get_params().c_side_name_translation
        if member.name in f_to_c_name:
            c_name = f_to_c_name[member.name]
        else:
            c_name = f_to_c_name.get(f"{fstruct.name}%{member.name}", member.name)

        return cls(
            is_component=True,
            f_name=member.name,
            c_name=c_name,
            type=type_.lower(),
            kind=member.kind or "",
            pointer_type=pointer_type,
            array=member.dimension.replace(" ", "").split(",") if member.dimension else [],
            init_value=str(member.default) if member.default else None,
            comment=member.comment,
            member=member,
        )

    @property
    def full_type(self):
        return FullType(self.type, len(self.array), self.pointer_type)

    @property
    def lbound(self) -> list[str]:
        if not self.array or self.array[0] == ":":
            return []

        return [dim.split(":")[0] if ":" in dim else "1" for dim in self.array]

    @property
    def ubound(self) -> list[str]:
        if not self.array or self.array[0] == ":":
            return []

        return [dim.split(":")[1] if ":" in dim else dim for dim in self.array]

    def get_dim1(self) -> tuple[str, str]:
        if self.ubound[0][-1] == "$":
            f_dim1 = self.ubound[0]
            c_dim1 = "Bmad::" + self.ubound[0][0:-1].upper()
            if self.lbound[0] != "1":
                logger.error('lbound not "1" with parameter upper bound!')
                sys.exit("STOPPING HERE")
            if self.ubound[0].lower() == "num_ele_attrib$":
                # NOTE: special case: this is an element attributes array, and we intend
                # to keep the array indices the same from C++/Fortran.
                return "num_ele_attrib$", "Bmad::NUM_ELE_ATTRIB+1"

        if not self.ubound[0].isnumeric():
            # NOTE: special case: n_pole_maxx->Bmad::N_POLE_MAXX
            return self.ubound[0], "Bmad::" + self.ubound[0].upper().rstrip("$")

        f_dim1 = str(1 + int(self.ubound[0]) - int(self.lbound[0]))
        c_dim1 = f_dim1
        return f_dim1, c_dim1

    @property
    def f_dims(self):
        if not self.array:
            return ()
        if len(self.array) == 1:
            return (self.f_dim1,)
        if len(self.array) == 2:
            return (self.f_dim1, self.dim2)
        if len(self.array) == 3:
            return (self.f_dim1, self.dim2, self.dim3)
        raise NotImplementedError(len(self.array))

    @property
    def c_dims(self):
        if not self.array:
            return ()
        if len(self.array) == 1:
            return (self.c_dim1,)
        if len(self.array) == 2:
            return (self.c_dim1, self.dim2)
        if len(self.array) == 3:
            return (self.c_dim1, self.dim2, self.dim3)
        raise NotImplementedError(len(self.array))

    @property
    def c_dim1(self) -> str:
        _, c_dim1 = self.get_dim1()
        return c_dim1

    @property
    def f_dim1(self) -> str:
        f_dim1, _ = self.get_dim1()
        return f_dim1

    @property
    def dim2(self) -> int:
        return 1 + int(self.ubound[1]) - int(self.lbound[1])

    @property
    def dim3(self) -> int:
        return 1 + int(self.ubound[2]) - int(self.lbound[2])

    def should_translate(self, struct_name: str) -> bool:
        params = get_params()
        return (
            self.kind not in params.component_no_translate_list
            and f"{struct_name}%{self.f_name}" not in params.component_no_translate_list
        )

    def _replace_transform_placeholders(self, transform: Transform) -> Transform:
        if self.type == "type":
            kind = self.kind
            if kind.lower().endswith("_struct"):
                kind = kind[: -len("_struct")]

            if not kind:
                raise RuntimeError("Kind is empty?")
            transform.replace_all("PROXYCLS", struct_to_proxy_class_name(kind))

        if self.pointer_type == NOT:
            # not a pointer/dynamically allocated type;
            # replace DIM1, DIM2, DIM3 here
            if len(self.array) >= 1 and self.array[0] not in {":", "0:"}:
                transform.replace_all("DIM1", str(self.c_dim1))
            if len(self.array) >= 2 and self.array[1] not in {":", "0:"}:
                transform.replace_all("DIM2", str(self.dim2))
            if len(self.array) >= 3 and self.array[2] not in {":", "0:"}:
                transform.replace_all("DIM3", str(self.dim3))
        return transform

    @property
    def transform(self):
        # self.is_dynamic_array
        transform = get_type_transform(self.full_type)
        return self._replace_transform_placeholders(transform)

    def original_repr(self) -> str:
        return f'["{self.type}({self.kind})", "{self.pointer_type}", "{self.f_name}", {self.array}, {self.lbound} {self.ubound} "{self.init_value}"]'


@dataclass
class CodegenStructure:
    f_name: str = ""  # Struct name on Fortran side
    short_name: str = ""  # Struct name without trailing '_struct'. Note: C++ name is 'CPP_<short_name>'
    cpp_class: str = ""  # C++ name.
    arg: list[Argument] = field(
        default_factory=list
    )  # List of structrure components + array bound dimensions.
    c_constructor_arg_list: str = ""
    c_constructor_body: str = ""  # Body of the C++ class_initializer
    c_extra_methods: str = ""  # Additional custom methods

    module: str = "unknown_module"
    parsed: ParsedStructure | None = None

    @property
    def container_alloc_name(self):
        return f"{self.f_name}_container_alloc"

    @property
    def args_to_convert(self):
        return [
            arg for arg in self.arg if f"{self.f_name}%{arg.f_name}" not in get_params().interface_ignore_list
        ]

    @property
    def recursive(self) -> bool:
        return any(
            arg.is_component
            and arg.type == "type"
            and arg.member is not None
            and arg.member.type_info.kind.lower() == self.f_name.lower()
            for arg in self.arg
        )

    def __str__(self) -> str:
        return f"[name: {self.short_name}, #arg: {len(self.arg)}]"
