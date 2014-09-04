#!/usr/bin/env python3

import sys
import re 
import json

# Python 2-3 compatability
try: input = raw_input
except NameError: pass

if len(sys.argv) <= 1:
    print("Usage: {} <file>...".format(sys.argv[0]))
    sys.exit(0)

with open("header_database.json") as f:
    headers = json.load(f)

for filename in sys.argv[1:]:
    with open(filename) as f:
        data = f.read()
    
    identifiers_used = set()
    for match in re.findall(r"std\s*::\s*([a-zA-Z_][a-zA-Z0-9_]*)", data):
        identifiers_used.add(match)
    
    headers_needed = set()
    for identifier_used in identifiers_used:
        if identifier_used not in headers:
            sys.stderr.write(identifier_used + ": ")
            headers[identifier_used] = input()
    
        headers_needed.add(headers[identifier_used])
    
    print("Includes for {}:".format(filename))
    for header_needed in sorted(headers_needed):
        print("#include <{}>".format(header_needed))
    print()

json.dump(headers, open("header_database.json", "w"), sort_keys=True, indent=4)

