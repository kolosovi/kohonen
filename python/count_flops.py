import sys

CODES = set([
    "r530210",
    "r538010",
    "r530410",
    "r531010",
    "r532010",
    "r534010",
    "r530810",
    "r530110"
])


def print_count(input):
    counts = 0

    for line in input:
        line = line.strip()

        try:
            count, code = line.split()
        except:
            continue
        else:
            if code.strip() in CODES:
                try:
                    count = int(count.strip())
                except:
                    continue
                else:
                    counts += count

    if counts > 0:
        print counts


if __name__ == "__main__":
    print_count(sys.stdin)
