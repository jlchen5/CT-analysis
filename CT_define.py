def determine_topology(contact1, contact2):
    c1_start, c1_end = sorted(contact1[1:3])
    c2_start, c2_end = sorted(contact2[1:3])
    if (c1_start < c2_start and c1_end < c2_end and c1_end < c2_start) or (c1_start < c2_start and c1_end < c2_end and c1_end == c2_start):
        return "S"
    elif (c1_start < c2_start and c1_end < c2_end ) or (c1_start < c2_start and c1_end < c2_end):
        return "P"
    elif (c1_start < c2_start and c2_start < c1_end ):
        return "X"
    else:
        return "N"

# path
base_input_path = '/Users/jialechen/Desktop/PhD/CT/Pang_2022_GenomeBiol_3D/data/hic_contact/'
base_output_path = '/Users/jialechen/Desktop/PhD/CT/Pang_2022_GenomeBiol_3D/results/topology_result/'

# read and analyze each file and save the results
for i in range(1, 25):
    input_file_path = f"{base_input_path}Cell_ID_{i:02d}.bed"  # 构建每个文件的路径
    output_file_path = f"{base_output_path}topology_result_Cell_ID_{i:02d}.txt"  # 输出文件的路径

    # read contacts
    contacts = []
    with open(input_file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            chr1, start1, end1, chr2, start2, end2 = parts[:6]
            contacts.append(((chr1, int(start1), int(end1)), (chr2, int(start2), int(end2))))

    # analyze contacts
    results = []
    for c1, c2 in contacts:
        topology = determine_topology(c1, c2)
        results.append(c1 + c2 + (topology,))

    # save results
    with open(output_file_path, 'w') as out_file:
        for result in results:
            line = f"{result[0]}\t{result[1]}\t{result[2]}\t{result[3]}\t{result[4]}\t{result[5]}\t{result[6]}\n"
            out_file.write(line)

    print(f"Results saved to {output_file_path}")
