import pysam
import numpy as np
import argparse


bin_default = 200


def snr_count(filepath, fragment_size, bin_size, p_low, p_high):
    counts = bin_coverage(filepath, fragment_size, bin_size)[1:]
    noise = get_quantile_from_freqs(p_low, counts)
    signal = get_quantile_from_freqs(p_high, counts)
    return (signal + 1) / (noise + 1)


def bin_coverage(filepath, fragment_size, bin_size):
    bins = []
    counts = [0 for _ in range(1000)]
    curr_ref = ""

    bf = pysam.AlignmentFile(filepath, "rb")

    for read in bf.fetch():
        if read.reference_name != curr_ref:
            bins_lens = [len(b) for b in bins]
            cs, ns = np.unique(bins_lens, return_counts=True)
            for i in range(len(cs)):
                if cs[i] >= len(counts):
                    for j in range(cs[i] - len(counts) + 1):
                        counts.append(0)
                counts[cs[i]] += ns[i]

            curr_ref = read.reference_name
            bins = [
                set()
                for _ in range(
                    bf.lengths[bf.references.index(read.reference_name)] // bin_size + 1
                )
            ]

        if read.is_reverse:
            pos = read.reference_end - fragment_size // 2
        else:
            pos = read.reference_start + fragment_size // 2

        try:
            bins[pos // bin_size].add(pos % bin_size)
        except IndexError:
            continue

    return counts


def get_quantile_from_freqs(q, freqs):
    idx = 0
    for i, f in enumerate(freqs):  # last nonzero element
        if f != 0:
            idx = i

    sum_ = sum(freqs[: idx + 1])
    i = -1

    cumsum = 0
    while cumsum / sum_ < q and i < idx:
        i += 1
        cumsum += freqs[i]

    return i


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Count signal-to-noise ratio in given bam file (index required)"
    )
    parser.add_argument(dest="filepath", type=str, help="file in .bam format")
    parser.add_argument("-d", dest="fragment_size", type=int, help="fragment size")
    parser.add_argument(
        "--bin_size",
        dest="bin_size",
        type=int,
        default=200,
        help=f"bin size for coverage counting (DEFAULT: {bin_default})",
    )
    parser.add_argument(
        "--plow",
        dest="p_low",
        type=float,
        default=0.1,
        help=f"Noise percentile (DEFAULT: 0.1)",
    )
    parser.add_argument(
        "--phigh",
        dest="p_high",
        type=float,
        default=0.9,
        help=f"Signal percentile (DEFAULT: 0.9)",
    )

    args = parser.parse_args()

    if args.fragment_size is None:
        raise Exception("fragment size is required")
    if args.p_low >= args.p_high:
        raise Exception("High percentile should be bigger that low one")

    print(snr_count(args.filepath, args.fragment_size, args.bin_size, args.p_low, args.p_high))
