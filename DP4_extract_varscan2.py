def varscan2_anno_DP4(dirpath,regexp):
    '''This function accepts the directory with all the vcf files to be opened. It also accepts the regex expression
    that will create the ID column based on the name of the file. E.g., 
    TCGA-123-DFB-23R2F_ANNOTATAED_VARSCAN2 can be collapsed to TCGA-123 as ID. The value returned is a dataframe with the
    chromosome, position of variation and the DP4 values '''
    #create the id_list which hold id values per file analyzed and the full path to invoke it for iterations.
    import os
    import re
    import pandas as pd
    DF = pd.DataFrame()
    filenames = os.listdir(dirpath)
    id_tcga = []
    full_path = []
    for filename in filenames:
        m = re.search(regexp, filename)
        id_tcga = m[1]
        full_path = (dirpath+filename)
        #open vcf file and pass it to a handler, then read it as a string
        with open(full_path, "rt") as f:
            vcf = f.read()

            #remove the INFO field. It is too concentrated with information to look at.
            def removeCSQ(line):
                    values = line.split()
                    INFO = values[7]
                    NEWINFO = []
                    for info in INFO.split(";"):
                            if not "CSQ=" in info:
                                    NEWINFO.append(info)
                    values[7] = ";".join(NEWINFO)
                    return "\t".join(values)

            #Detect the field I am interested and group them in endlist
            line = vcf.split()
            line[11]

            n = 11
            endlist = [[] for _ in range(n)]
            for index, item in enumerate(line):
                endlist[index % n].append(item)


            # list splice to create a dp4 list
            endlist[-1]
            n = -1
            dp4_abcd = []
            for i in endlist[-1]:
                n = n + 1
                dp4 = endlist[-1][n].split("%:")
                dp4 = dp4[1]
                dp4_abcd.append(dp4)
            # split the dp4 string to a list    
            n = -1
            for i in dp4_abcd:
                n = n + 1
                dp4_abcd[n] = dp4_abcd[n].split(",")

            # match 1-1 elements from multiple lists to create a list of lists     
            comb_lists = [[x] + [y] +[z] for x, y, z in zip(endlist[0], endlist[1], dp4_abcd)]
            #print(comb_lists) 
            #convert list of lists to pandas dataframe
            columns = ["chrom","pos","abcd"]
            df = pd.DataFrame(comb_lists, columns = columns)
            #create a,b,c,d columns for ref fwd, ref rev, alt fwd, alt rev
            n = -1
            df["a"] = "NA"
            df["b"] = "NA"
            df["c"] = "NA"
            df["d"] = "NA"
            for i in df.iterrows():
                n = n + 1
                df["a"][n] = df["abcd"][n][0]
                df["b"][n] = df["abcd"][n][1]
                df["c"][n] = df["abcd"][n][2]
                df["d"][n] = df["abcd"][n][3]
            #remove the combined abcd
            df.drop('abcd', inplace=True, axis=1)
            #add a filename idendifier
            df["ID"] = id_tcga
            DF = DF.append(df)
            DF = DF.astype({"a": int, "b": int, "c": int, "d": "int"})
    return DF
