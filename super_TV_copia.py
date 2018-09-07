# NERIS GARCIA GONZALEZ
import sys
import getopt
from collections import Counter
import shutil
import os
import csv
import operator


# Filtering Blast results keeping only the best hit and those hits with
# a Bit Score higher than 90 % of the best hit Bit Score.
def filter_by_bitscore(ablast_file, sample_name):
    blast_file = open(ablast_file, "r")
    out_file = open(sample_name + "_filtered.txt", "w")
    # Iniciate variables reading first line
    line1 = blast_file.readline()
    lline1 = line1.split()
    # Iniciate first gene and its bitscrore
    gene = str(lline1[0])
    bitscore = float(lline1[11])
    # We always keep the first line, the besthit
    out_file.write(str(line1))
    # Calculate bitscore threshold
    bitscore10 = bitscore - (bitscore * 10 / 100)
    # Start reading the rest of the file
    for line in blast_file:
        lline = line.split()
        # Check if we still looking at the same gene
        if str(lline[0]) == gene:
            if float(lline[11]) >= bitscore10:
                # If a hit has a Bit Score higher than 90 % of the best hit
                # Bit Score it is wrote in the output
                out_file.write(str(line))
        # If we are not looking at the same gene, write ouput and
        # re-start variables
        else:
            # We always keep the first line, the besthit.
            out_file.write(str(line))
            gene = str(lline[0])
            bitscore = float(lline[11])
            bitscore10 = bitscore - (bitscore * 10 / 100)


# Add to filtered Blast results, the taxonomy of the hit
def blast2taxonomy(blast_filtered_file, db_file, sample_name):
    blast_file = open(blast_filtered_file, "r")
    out = open(sample_name + "_filtered_tax", "w")
    # Creation of a dictionary with the taxonomy database. The key will be
    # the key of the taxonomy and the values the taxonomic range.
    # BctDB = open("BacteriaDB_2016_08_08.taxonomy", "r")
    BctDB = open(db_file, "r")
    dict_db = {}
    for line in BctDB:
        if not line.startswith("Plasmids"):
            num = line.split("Bacteria")[0]
            tax = line.split("Bacteria")[1]
            dict_db[num[:-1]] = "Bacteria\t" + tax
    # Creation of a dictionary with Blast results. As a key, the key of the
    # taxonomy and as value, the blast results (all the line)
    dict_blast = {}
    for line in blast_file:
        if line.endswith("\n"):
            line = line[:-1]
        numero = line.split(".")[0].split()[1]
        dict_blast[line] = numero
    # Cross of two dictoraries to obtain Bast hits and its taxonomy
    for k, v in dict_blast.iteritems():
        # Write results in an external file
        out.write(k + "\t" + dict_db[v])
        print k + "\t" + dict_db[v]


# Apply Taxonomic Vote algorithin to filtered Blasts results with taxonomy
def taxonomic_vote(blast_file, sample_name):
    output = open(sample_name + ".TV.out", "w")
    # Diccionario donde guardaremos el hit consenso para cada gen
    output.write(
        "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\t"
        "sstart\tsend\tevalue\tbitscore\tkingdom\tphylum\tclass\torder\t"
        "family\tgenus\tspecies\tTV\n")
    dhits = {}
    # Iniciate variables reading first line
    llinea1 = blast_file.readline().split()
    gene = str(llinea1[0])
    besthit = llinea1
    # All the hits of a gene will be keept in a matrix. All lines for
    # a gene will be rows of the matrix.
    matrix = []
    matrix.append(llinea1)
    # Start reading the file
    for line in blast_file:
        # Adding all the lineas for a same gene to its matrix
        if str(line.split()[0]) == str(gene):
            matrix.append(line.split())
        # When matrix is complete, it is processed.
        else:
            # TV_type is the Taxomoic Vote class for that gene
            # Only one Best Hit
            if len(matrix) <= 1:
                TV_type = "I\tSp."
            # More than one hit
            else:
                # Taxonomic Vote algorithym is applied
                # for each taxon level, number of occcurrences are counted
                # Variable y will mark the taxon level
                y = 0
                for taxon_level in ["Sp.", "Genus", "Family", "Order",
                                    "Class", "Phylum", "Domain"]:
                    matrix_taxon_level = []
                    # If we are at the species level, we have to save from
                    # element 18 of the list, until the last element
                    if y == 0:
                        bhit = besthit[18:]
                        for hit in matrix:
                            # For each hit, the corresponding taxon level,
                            # is added to a matrix
                            matrix_taxon_level.append(str(hit[18:]))
                    # If we are not at the species level, we only keep
                    # a single element
                    else:
                        bhit = besthit[18 - y]
                        for hit in matrix:
                            matrix_taxon_level.append(str(hit[18 - y]))
                    # Creation of a dictionary with the number of
                    # appeareances for each taxon
                    dTaxonCount = Counter(matrix_taxon_level)
                    # Counting the number of appeareances of the Best Hit
                    N_bhit = dTaxonCount.get(str(bhit))
                    # With the number of appearances of each taxon we can
                    # know the Taxonomic Vote type of the gene
                    # If the best hit es the major represented:
                    if float(N_bhit) >= 0.5 * len(matrix_taxon_level):
                        if N_bhit == 0.5 * len(matrix_taxon_level):
                            lista_matrix_nobhit_set = set()
                            for mtl in matrix_taxon_level:
                                lista_matrix_nobhit_set.add(mtl)
                            lista_matrix_nobhit_set.remove(str(bhit))
                            hit2 = (lista_matrix_nobhit_set.pop())
                            N_hit2 = dTaxonCount.get(str(hit2))

                            # If N_besthit == N_hit2 it could be a tie
                            if N_hit2 == 0.5 * len(matrix_taxon_level):
                                tie = set()
                                for elemento in matrix:
                                    tie.add(elemento[17])

                                if len(tie) == 1:
                                    TV_type = "II\t" + str(taxon_level)
                                    break

                                else:
                                    x = 0
                                    for taxa in ["Family", "Order",
                                                 "Class", "Phylum",
                                                 "Domain"]:
                                        tie = set()
                                        for elemento in matrix:
                                            tie.add(elemento[16 - x])
                                        if len(tie) == 1:
                                            if x == 4:
                                                TV_type = "V\t" + str(taxa) + "\t" + str(hit2).replace("[","").replace("]", "").replace(",","").replace("\'", "")
                                            else:
                                                TV_type = "III\t" + str(taxa)
                                            break
                                        x += 1
                                    break
                        if y == 0:
                            TV_type = "II\tSp."
                            break
                        elif y == 1:
                            TV_type = "II\tGenus"
                            break
                        else:
                            TV_type = "III\t" + str(taxon_level)
                            break
                    # If the besthit is not the major taxon represented
                    # we have to obtain the major one
                    if N_bhit < 0.5 * len(matrix_taxon_level):
                        if max(dTaxonCount.values()) > 0.5 * len(matrix_taxon_level):
                            major_hit = [key for key, value in dTaxonCount.iteritems()
                                 if value == max(dTaxonCount.values())][0]
                            TV_type = "IV\t" + str(taxon_level) + "\t" + str(major_hit).replace("[", "").replace("]", "").replace(",", "").replace("\'", "")
                            break
                    y += 1
            # Write output for this gene
            hit_out = ""
            for bh in besthit[:18]:
                hit_out += str(bh) + "\t"
            output.write(hit_out + str(besthit[18:]).replace("[", "").replace("]","").replace(",", "").replace("\'", "") + "\t" + TV_type + "\n")
            # Re-start variables for next gene
            gene = str(line.split()[0])
            besthit = line.split()
            matrix = []
            matrix.append(line.split())
            TV_type = ""


# To know how the genes are located in the contigs ans its taxonomy.
def genes_per_contig(gff):
    # Create a dictionary were the key will be the contig name and the
    # value the genes inside the contig
    dContigGenes = {}
    # open the annotation file
    gff_file = open(gff, "r")
    for line in gff_file:
        # Skipping coments and fasta sequences
        if line.startswith(">"):
            break
        if not line.startswith("##"):
            # Reading files with CDS with the annotation of each gene
            if "CDS" in line:
                lline = line.split()
                if lline[0] not in dContigGenes:
                    lgenes = []
                    lgenes.append(lline)
                    dContigGenes[lline[0]] = lgenes
                else:
                    lgenes.append(lline)
                    dContigGenes[lline[0]] = lgenes
    return dContigGenes


# Convert Blast results to a dictionary
def blast2dict(TV_file):
    dTV = {}
    TV = open(TV_file, "r")
    for line in TV:
        dTV[str(line.split("\t")[0])] = (line.split("\t")[1:])
    TV.close()
    return dTV


def genes_unicos(dBlast, gff_file, sample_name):
    dict = {}
    HP = 0
    list = []
    gff = open(gff_file, "r")
    out = open(sample_name + "_gene_info.tab", "w")
    for line in gff:
        # primero miramos los CDS (genes identificados con UbiprotKB, HAMAP o Hypothetical Proteins)
        if "CDS" in line:

            if "UniProtKB" in line:
                # Miramos si esta en el dict
                UniprotGeneID = "UniProtKB:" + str(
                    line.split("UniProtKB")[1].split(";")[0][1:])
                if UniprotGeneID not in dict:
                    dict[UniprotGeneID] = 1
                else:
                    dict[UniprotGeneID] = int(dict[UniprotGeneID]) + 1

                    # Creamos la lista con la info

                nlist = []
                # Uniprot
                nlist.append(UniprotGeneID)
                # Producto
                nlist.append(line.split("product=")[1][:-1])
                ProkkageneID = line.split(";")[0].split("=")[1]
                nlist.append(ProkkageneID)
                # contig
                nlist.append(line.split("\t")[0])
                nlist.append(line.split("\t")[3])
                nlist.append(line.split("\t")[4])

                # datos blast desde dBlast
                if ProkkageneID in dBlast:
                    nlist.append(
                        (dBlast[line.split(";")[0].split("=")[1]])[1])
                    nlist.append(
                        (dBlast[line.split(";")[0].split("=")[1]])[2])
                    nlist.append(
                        (dBlast[line.split(";")[0].split("=")[1]])[10:])

                list.append(nlist)



            elif "HAMAP" in line:
                hamap = "HAMAP" + str(line.split("HAMAP")[1].split(";")[0])

                if hamap not in dict:
                    dict[hamap] = 1
                else:
                    dict[hamap] = int(dict[hamap]) + 1

                # list.append()#gen
                # print line.split("UniProtKB")[1].split(";")[0][1:]
                nlist = []
                nlist.append(hamap)
                nlist.append(line.split("product=")[1][:-1])
                ProkkageneID = line.split(";")[0].split("=")[1]
                nlist.append(ProkkageneID)
                # contig
                nlist.append(line.split("\t")[0])
                nlist.append(line.split("\t")[3])
                nlist.append(line.split("\t")[4])

                # datos blast desde dBlast
                if ProkkageneID in dBlast:
                    nlist.append(dBlast[ProkkageneID][1])
                    nlist.append(dBlast[ProkkageneID][2])
                    nlist.append(dBlast[ProkkageneID][10:])

                list.append(nlist)


            elif "hypothetical protein" in line:
                HP += 1
                HP_Name = "Hypothetical protein " + str(HP)
                dict[HP_Name] = "-"
                nlist = []
                # Do not have uniprot id
                nlist.append(HP_Name)
                # Do not have Product
                nlist.append("-")
                # Prokka id
                ProkkageneID = line.split(";")[0].split("=")[1]
                nlist.append(ProkkageneID)
                # Contig information
                nlist.append(line.split("\t")[0])  # Contig
                nlist.append(line.split("\t")[3])  # Contig start
                nlist.append(line.split("\t")[4])  # Contig end
                # Add Blast information
                if ProkkageneID in dBlast:
                    nlist.append(dBlast[ProkkageneID][1])
                    nlist.append(dBlast[ProkkageneID][2])
                    nlist.append(dBlast[ProkkageneID][10:])
                list.append(nlist)

                # Si no es un CDS es porque es rRNA o tRNA
        elif "barrnap" in line:  # rRNA
            barrnap = line.split("product=")[1][:-1]

            if barrnap not in dict:
                dict[barrnap] = 1
            else:
                dict[barrnap] = dict[barrnap] + 1

            nlist = []
            nlist.append(barrnap)
            nlist.append("rRNA")
            ProkkageneID = line.split(";")[0].split("=")[1]
            nlist.append(ProkkageneID)
            # contig
            nlist.append(line.split("\t")[0])
            nlist.append(line.split("\t")[3])
            nlist.append(line.split("\t")[4])

            # datos blast desde dBlast
            if ProkkageneID in dBlast:
                nlist.append(dBlast[ProkkageneID][1])
                nlist.append(dBlast[ProkkageneID][2])
                nlist.append(dBlast[ProkkageneID][10:])

            list.append(nlist)

        elif "Aragorn" in line:  # tRNA
            Aragorn = line.split("product=")[1][:-1]
            if Aragorn not in dict:
                dict[Aragorn] = 1
            else:
                dict[Aragorn] = dict[Aragorn] + 1

            nlist = []
            nlist.append(Aragorn)
            if "tmRNA" in line:  # Puede haber tmRNA
                nlist.append("tmRNA")
            else:
                nlist.append("tRNA")
            ProkkageneID = line.split(";")[0].split("=")[1]
            nlist.append(ProkkageneID)
            # contig
            nlist.append(line.split("\t")[0])
            nlist.append(line.split("\t")[3])
            nlist.append(line.split("\t")[4])

            # datos blast desde dBlast
            if ProkkageneID in dBlast:
                nlist.append(dBlast[ProkkageneID][1])
                nlist.append(dBlast[ProkkageneID][2])
                nlist.append(dBlast[ProkkageneID][10:])

            list.append(nlist)

    # Para imprimir lo anterior en pantalla
    out.write(
        "DataBase:ID\tRepetitions\tProduct\tPROKKA gene ID\tContig\tend"
        "\tstart\tpident\tlength\tbitscore\tkingdom\tphylum\tclass\torder"
        "\tfamily\tgenus\tspecies\tTV\n")
    for k, v in dict.iteritems():
        for lista in list:
            if str(k) in lista:
                txt = str(k) + "\t" + str(v) + "\t" + str(
                    lista[1:]).replace("[", "").replace("\\n", "").replace(
                    "]", "").replace("\'", "").replace(",", "\t") + "\n"
                out.write(txt)

    gff.close()
    # print list


def main(argv):
    # Argumentos...........................................................
    blast_file = ''
    gff_file = ''
    db_file = ''
    try:
        opts, args = getopt.getopt(argv, "hd:b:a:",
                                   ["DataBase=", "blast_file=",
                                    "gff_file="])
    except getopt.GetoptError:
        print "\nsuper_TV.py -d <Database Taxonomy> -b <blast_file> " \
              "-a <annotation_file>\n\t-d\tDatabase: *.taxonomy\n\t" \
              "-b\tBlast output \n\t-a\tAnnotation file format: gff\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print "\nsuper_TV.py -d <Database Taxonomy> -b <blast_file> " \
                  "-a <annotation_file>\n\t-d\tDatabase: *.taxonomy\n\t" \
                  "-b\tBlast output \n\t-a\tAnnotation file format: gff\n"
            sys.exit()
        elif opt in ("-d", "--dbfile"):
            db_file = arg
        elif opt in ("-b", "--ifile"):
            blast_file = arg
        elif opt in ("-a", "--gffile"):
            gff_file = arg
    # .....................................................................
    sample_name = \
        blast_file.split("/")[len(blast_file.split("/")) - 1].split("_")[0]
    filter_by_bitscore(blast_file, sample_name)
    blast_filered_file = sample_name + "_filtered.txt"
    blast2taxonomy(blast_filered_file, db_file, sample_name)

    # Para ordenar el fichero necesitamos una copia
    blast_filered_tax_file_name = sample_name + "_filtered_tax"

    blast_filered_tax_file2 = sample_name + "_filtered_tax2"
    shutil.copy(blast_filered_tax_file_name, blast_filered_tax_file2)

    blast_filered_tax_file = open(blast_filered_tax_file_name, "w")

    # Ordenamos el fochero
    reader = csv.reader(open(blast_filered_tax_file2), delimiter="\t")
    for line in sorted(reader, key=operator.itemgetter(0)):
        cat = ""
        for elemento in line:
            cat += str(elemento) + "\t"
        cat += "\n"
        blast_filered_tax_file.write(cat)

    blast_filered_tax_file.close()
    os.remove(blast_filered_tax_file2)

    #   ELIMINAR COLUMNA DE MAS
    # abrimos el fichero con los resultados de  blast
    shutil.copy(blast_filered_tax_file_name, "blast_file_addfinalgene")
    blast_file = open("blast_file_addfinalgene", "a")
    blast_file.write("Hyh2FghGH254kjuiIF5FJKJsgf31dFJ")
    blast_file.close()
    blast_file = open("blast_file_addfinalgene", "r")

























    dhits = taxonomic_vote(blast_file, sample_name)
    os.remove("blast_file_addfinalgene")

    dctg = genes_per_contig(gff_file)
    TV_file_name = sample_name + ".TV.out"
    dBlast = blast2dict(TV_file_name)
    genes_unicos(dBlast, gff_file, sample_name)


if __name__ == '__main__':
    main(sys.argv[1:])
