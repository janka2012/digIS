import csv
import logging
import os


def read_csv(filename, delimiter=';'):
    out_list = []
    if os.path.exists(filename):
        with open(filename, newline='') as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            for row in reader:
                out_list.append(row)
    else:
        logging.error("Filename {} does not exist.".format(filename))
    return out_list


def write_csv(rows, filename, header, delimiter=","):
    if len(rows) > 0 and len(header) != len(rows[0]):
        raise ValueError("Number of elements in header and in row to write is not same.")\

    with open(filename, 'w+', newline='') as f:
        writer = csv.writer(f, delimiter=delimiter)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
