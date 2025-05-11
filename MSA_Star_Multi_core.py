import time
import numpy as np
import multiprocessing as mp

from PSA_Kband import PSA_AGP_Kband
from score import spscore
from FASTA import readfasta


def find_center_sequence_sub(task_name: str, parameters: list):
    strs = parameters[0]
    task_id = parameters[1]
    total_tasks = parameters[2]
    length = len(strs)

    chunk_size = (length * (length - 1) // 2) // total_tasks + 1
    results = [0] * chunk_size

    for i in range(length):
        for j in range(i + 1, length):
            idx = calculate_index(length, i, j)
            if idx % total_tasks == task_id:
                score, _ = PSA_AGP_Kband(A=strs[i], B=strs[j], get_score=1)
                results[idx // total_tasks] = score

    return results


def calculate_index(length: int, i: int, j: int):
    return ((2 * length - i - 1) * i // 2) + j - i


def add_gaps(mark: list, sequence: str):
    result = ""
    for i in range(len(mark)):
        result += "-" * mark[i]
        if i < len(mark) - 1:
            result += sequence[i]
    return result


def find_center_sequence(strs: list, cores: int):
    length = len(strs)
    distance_matrix = [[0] * length for _ in range(length)]
    tasks = {}
    for i in range(cores):
        tasks[f"task{i}"] = (strs, i, cores)
    pool = mp.Pool(cores)
    results = [pool.apply_async(find_center_sequence_sub, args=(name, param)) for name, param in tasks.items()]
    results = [p.get() for p in results]

    for i in range(length):
        for j in range(i + 1, length):
            idx = calculate_index(length, i, j)
            distance_matrix[i][j] = distance_matrix[j][i] = results[idx % cores][idx // cores]

    return np.argmax(np.sum(distance_matrix, axis=0))


def pairwise_align(strs: list, center_idx: int):
    aligned_strs = []
    for i in range(len(strs)):
        if i != center_idx:
            _, aligned_center, aligned_other = PSA_AGP_Kband(strs[center_idx], strs[i])
            aligned_strs.append([aligned_center, aligned_other])
    return aligned_strs


def calculate_gap_locations(aligned_strs: list, gap_marks: list, center_idx: int):
    for aligned_str in aligned_strs:
        i = 0
        gap_count = 0
        for char in aligned_str[0]:
            if char == '-':
                gap_count += 1
            else:
                gap_marks[i] = max(gap_marks[i], gap_count)
                gap_count = 0
                i += 1
        gap_marks[i] = max(gap_marks[i], gap_count)
    return gap_marks


def insert_gaps_in_sequences(aligned_strs: list, gap_marks: list, strs: list, center_idx: int):
    aligned_sequences = [""] * len(strs)
    aligned_sequences[center_idx] = add_gaps(gap_marks, strs[center_idx])
    idx = 0
    for aligned_str in aligned_strs:
        mark = [0] * (len(aligned_str[0]) + 1)
        total = 0
        pi = 0
        pj = 0
        for char in aligned_str[0]:
            if char == '-':
                total += 1
            else:
                mark[pi] = gap_marks[pj] - total
                pi += 1
                pj += 1
                while total != 0:
                    pi += 1
                    total -= 1
        mark[pi] = gap_marks[pj] - total
        if idx >= center_idx:
            aligned_sequences[idx + 1] = add_gaps(mark, aligned_str[1])
        else:
            aligned_sequences[idx] = add_gaps(mark, aligned_str[1])
        idx += 1
    return aligned_sequences


def run_multicore_msa(strs, cores=-1):
    start_time = time.time()

    if cores == -1:
        cores = mp.cpu_count()

    if cores is None:
        raise ValueError("Unable to detect number of cores. Please specify the 'cores' parameter.")

    print("-----------RUN-----------")
    print("Using cores:", cores)
    print("Loading center sequence...", end="")

    center_idx = find_center_sequence(strs, cores)

    print("Center sequence loaded.")
    print("Center sequence:", ''.join(strs[center_idx]))

    aligned_strs = pairwise_align(strs, center_idx)

    gap_marks = [0] * (len(strs[center_idx]) + 1)
    gap_marks = calculate_gap_locations(aligned_strs, gap_marks, center_idx)
    aligned_strs = insert_gaps_in_sequences(aligned_strs, gap_marks, strs, center_idx)

    score_value = spscore(aligned_strs)
    end_time = time.time()

    print("Execution time: %.2f seconds" % (end_time - start_time))
    print("SP score:", score_value)

    for aligned_str in aligned_strs:
        print(' '.join(aligned_str))

    print("-----------END-----------")

    return score_value, aligned_strs
