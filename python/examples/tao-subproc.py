#!/usr/bin/env python3
from __future__ import annotations

from pytao import SubprocessTao
from tao_subproc_helper import get_info_example


def test():
    with SubprocessTao(
        init_file="$ACC_ROOT_DIR/bmad-doc/tao_examples/optics_matching/tao.init",
        noplot=True,
    ) as tao:
        print("Tao is:", tao)
        print("\n".join(tao.cmd("show ele 1")[:10]))

        data = tao.subprocess_call(get_info_example)
        print()
        print("Retrieved info from subprocess with pybmad:")
        print(f"data = {data}")


if __name__ == "__main__":
    test()
