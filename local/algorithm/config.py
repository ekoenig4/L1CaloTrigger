from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-debug',action='store_true',default=False)

config = parser.parse_args()
