#!/usr/bin/env python
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

from __future__ import print_function
import sys
import os

import ananse.network


def network(args):

    if not os.path.exists(args.features):
        print("File %s does not exist!" % args.features)
        sys.exit(1)

    b = ananse.network.Network()
    b.run_network(args.features, args.outfile)
