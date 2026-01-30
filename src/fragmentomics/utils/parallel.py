"""
Parallel processing utilities.
"""

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from typing import Callable, Iterable, TypeVar, Optional, Sequence
import multiprocessing as mp

T = TypeVar("T")
R = TypeVar("R")


def parallel_map(
    func: Callable[[T], R],
    items: Sequence[T],
    n_workers: Optional[int] = None,
    use_threads: bool = False,
    show_progress: bool = True,
) -> list[R]:
    """
    Apply a function to items in parallel.
    
    Parameters
    ----------
    func : Callable
        Function to apply to each item
    items : Sequence
        Items to process
    n_workers : int, optional
        Number of workers (default: CPU count)
    use_threads : bool, default False
        Use threads instead of processes
    show_progress : bool, default True
        Show progress bar
        
    Returns
    -------
    list
        Results in same order as input items
    """
    if n_workers is None:
        n_workers = mp.cpu_count()
    
    n_workers = min(n_workers, len(items))
    
    if n_workers <= 1:
        # Serial execution
        if show_progress:
            try:
                from tqdm import tqdm
                return [func(item) for item in tqdm(items, desc="Processing")]
            except ImportError:
                pass
        return [func(item) for item in items]
    
    Executor = ThreadPoolExecutor if use_threads else ProcessPoolExecutor
    
    with Executor(max_workers=n_workers) as executor:
        if show_progress:
            try:
                from tqdm import tqdm
                results = list(tqdm(
                    executor.map(func, items),
                    total=len(items),
                    desc="Processing",
                ))
            except ImportError:
                results = list(executor.map(func, items))
        else:
            results = list(executor.map(func, items))
    
    return results
