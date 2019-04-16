def write_gff(rows, filename, header):
    if len(rows) > 0 and len(header) != len(rows[0]):
        raise ValueError("Number of elements in header and in row to write is not same.")\

    sid_idx = header.index('sid')
    start_idx = header.index('sstart')
    end_idx = header.index('send')
    strand_idx = header.index('strand')
    score_idx = header.index('acc')
    attr_idx = [x for x in range(len(header)) if x not in [sid_idx, start_idx, end_idx, strand_idx, score_idx]]

    with open(filename, 'w', newline='') as f:
        f.write("##gff-version 3\n")
        for row in rows:
            sid = row[sid_idx]
            start = row[start_idx]
            end = row[end_idx]
            strand = row[strand_idx]
            score = row[score_idx]
            attr_list = []
            for idx in attr_idx:
                attr_list.append(header[idx] + "=" + str(row[idx]))
            attributes = ";".join(attr_list)

            f.write("{}\tdigIS\ttransposable_element\t{}\t{}\t{}\t{}\t.\t{}\n".format(sid,start, end, score, strand, attributes))
        f.close()

# https://www.ensembl.org/info/website/upload/gff3.html
# bedtools getfasta -fi input.fasta -bed input.gff -fo output.fasta
# bedtools flank -i <BED/GFF/VCF> -g <GENOME> [-b or (-l and -r)]