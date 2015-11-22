#!/usr/bin/env python3

import sys
import re
import os
import readline # NOQA

# Python 2-3 compatability.
try:
    input = raw_input
except NameError:
    pass

try:
    os_replace = os.replace
except AttributeError:
    def os_replace(src, dst):
        if os.name == "nt":
            os.unlink(dst)
        os.rename(src, dst)


IWYU_PATH = os.path.dirname(os.path.realpath(__file__))
VALID_WILDCARD = re.compile(r"(?:[*a-zA-Z_][*a-zA-Z0-9_]*\s*::\s*)*[*a-zA-Z_][*a-zA-Z0-9_]*")
IDENTIFIERS = re.compile(r"\b(?:[a-zA-Z_][a-zA-Z0-9_]*\s*::\s*)*[a-zA-Z_][a-zA-Z0-9_]*")
STRIP_WS = re.compile(r"\s*")


def validate_wildcard(wildcard, linenr=0):
    wildcard = wildcard.strip()
    if not re.match(VALID_WILDCARD, wildcard):
        raise RuntimeError("invalid wildcard on line {}".format(linenr + 1))
    return wildcard.replace("*", ".*")


def make_wildcards_regex(wildcards):
    return re.compile("^{}$".format("|".join("(?:(){})".format(w[1]) for w in wildcards)))


def load_database():
    headers = {}
    wildcards = []
    with open(os.path.join(IWYU_PATH, "header_database.txt")) as db:
        for linenr, line in enumerate(db):
            line = line.strip()
            if line.startswith('?'):
                wildcards.append(["?", validate_wildcard(line[1:], linenr)])
            elif line.startswith('!'):
                wildcards.append(["!", validate_wildcard(line[1:], linenr)])
            elif '=' in line:
                lhs, rhs = line.split('=')
                if '*' in lhs:
                    wildcards.append(["=", validate_wildcard(lhs, linenr), rhs.strip()])
                else:
                    headers[lhs.strip()] = rhs.strip()
            elif line:
                raise RuntimeError("database format error on line {}".format(linenr=1))

    wildcards.sort(key=lambda w: w[0] == "?")

    return headers, wildcards


def store_database(headers, wildcards):
    tmp = os.path.join(IWYU_PATH, "header_database.txt~")
    dst = os.path.join(IWYU_PATH, "header_database.txt")

    with open(tmp, "w") as db:
        wildcards.sort(key=lambda w: w[0] != "?")
        for w in wildcards:
            if w[0] in "?!":
                db.write("{} {}\n".format(w[0], w[1].replace('.*', '*')))
            else:
                db.write("{} = {}\n".format(w[1].replace('.*', '*'), w[2]))

        for k, v in sorted(headers.items(), key=lambda kv: kv[::-1]):
            db.write("{} = {}\n".format(k, headers[k]))

    os_replace(tmp, dst)


if len(sys.argv) <= 1:
    print("Usage: {} <file>...".format(sys.argv[0]))
    sys.exit(0)

headers, wildcards = load_database()
wildcards_regex = make_wildcards_regex(wildcards)

for filename in sys.argv[1:]:
    with open(filename) as f:
        data = f.read()

    headers_needed = set()
    for identifier in re.findall(IDENTIFIERS, data):
        identifier = re.sub(STRIP_WS, "", identifier)
        if identifier in headers:
            headers_needed.add(headers[identifier])
            continue

        m = re.match(wildcards_regex, identifier)
        if m:
            wildcard = wildcards[m.lastindex - 1]
            if wildcard[0] == "!":
                continue
            elif wildcard[0] == "=":
                headers_needed.add(wildcard[2])
            else:
                sys.stderr.write(identifier + ": ")
                header = input().strip()
                if header == "!":
                    wildcards.append(["!", header])
                    wildcards_regex = make_wildcards_regex(wildcards)
                else:
                    headers[identifier] = header
                    headers_needed.add(header)

                store_database(headers, wildcards)

    print("Includes for {}:".format(filename))
    # Key is a hack to make " sort after <.
    for header_needed in sorted(headers_needed, key=lambda h: h.replace('"', "~")):
        print("#include {}".format(header_needed))
    sys.stdout.write("\n")
    sys.stdout.flush()

store_database(headers, wildcards)
