"""Tasks that misbehave on purpose, for tests/test_batch.py."""
import ctypes
import os
import time


def behave(path, **kwargs):
    mode = os.path.basename(path)
    if mode.startswith("ok"):
        return {"mode": mode, **kwargs}
    if mode.startswith("raise"):
        raise ValueError("deliberate failure")
    if mode.startswith("abort"):
        os.abort()
    if mode.startswith("segfault"):
        ctypes.string_at(0)
    if mode.startswith("hang"):
        time.sleep(600)
    return None
