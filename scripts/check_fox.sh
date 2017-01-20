#!/bin/sh

python -c "
import math
q = math.sqrt($1)
if q * q == $1 and $2 % q == 0: print True
else: print False
"
