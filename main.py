import sys
import numpy as np
import strawC
from scipy import sparse


def get_chromosomes(filepath: str):
    chrom_dot_sizes = strawC.getChromosomes(filepath)
    chromosomes = []
    for chromosome in chrom_dot_sizes:
        chrom = chromosome.name
        if chrom.lower() == 'all':
            continue
        chromosomes.append(chrom)
    return chromosomes


def get_sparse_matrix(filename: str, chrom1: str, chrom2: str, res: int):
    records = strawC.strawC("observed", "NONE", filename, chrom1, chrom2, "BP", res)
    row_indices, col_indices, data = list(), list(), list()
    for record in records:
        row_indices.append(record.binX)
        col_indices.append(record.binY)
        data.append(record.counts)
    row_indices = (np.asarray(row_indices)) / res
    col_indices = (np.asarray(col_indices)) / res
    r_width = int(np.max(row_indices)) + 1
    c_width = int(np.max(col_indices)) + 1
    data = sparse.coo_matrix((data, (row_indices.astype(int), col_indices.astype(int))),
                             shape=(r_width, c_width)).toarray()
    data[np.isnan(data)] = 0
    data[np.isinf(data)] = 0
    return data


def get_difference_between_matrices(file1: str, file2: str, chrom1: str, chrom2: str, res: int):
    matrix1 = get_sparse_matrix(file1, chrom1, chrom2, res)
    matrix2 = get_sparse_matrix(file2, chrom1, chrom2, res)
    return abs(matrix1 - matrix2).sum()


def get_norm_vector(filename: str, chrom: str, norm: str, res: int):
    try:
        footer = strawC.getNormExpVectors(filename, chrom, chrom, "observed", norm, "BP", res)
        if len(footer.c1Norm) > 1:
            data = np.asarray(footer.c1Norm)
            data[np.isnan(data)] = 0
            data[np.isinf(data)] = 0
            return data
    except:
        pass
    return 0



def get_difference_between_vector(file1, file2, chrom, norm, res):
    vector1 = get_norm_vector(file1, chrom, norm, res)
    vector2 = get_norm_vector(file2, chrom, norm, res)
    diff = np.abs(vector1 - vector2)
    return np.sum(diff), np.mean(diff), np.max(diff)


class CompareFiles:
    def __init__(self, file1: str, file2: str, norms: list):
        self.file1 = file1
        self.file2 = file2
        self.resolutions = strawC.getResolutions(file1)
        self.chromosomes = get_chromosomes(file1)
        self.norms = norms

    def check_every_chromosome_pair(self):
        n = len(self.chromosomes)
        total_error = 0
        for res in self.resolutions:
            res_error = 0
            for i in range(n):
                for j in range(i, n):
                    try:
                        res_error += get_difference_between_matrices(self.file1, self.file2, self.chromosomes[i],
                                                                self.chromosomes[j], res)

                        #print("CHR", self.chromosomes[i], "CHR", self.chromosomes[j],
                        #      "RES", res, "ERROR", error)
                    except:
                        pass
            print("Total error between entries for resolution ", res, ":", res_error)
            total_error += res_error
        print("Total error for matrices across all resolutions ", total_error)

    def check_all_normalizations(self):
        total_error, total_mu_error, total_max_error = 0, [], []
        for norm in self.norms:
            norm_error, norm_mu_error, norm_max_error = 0, [], []
            for chrom in self.chromosomes:
                for res in self.resolutions:
                    err, mu_err, max_err = get_difference_between_vector(self.file1, self.file2, chrom, norm, res)
                    norm_error += err
                    norm_mu_error.append(mu_err)
                    norm_max_error.append(max_err)
            #print("NORM", norm, "ERROR", norm_error, "MEAN ERROR", np.mean(norm_mu_error))
            total_error += norm_error
            total_mu_error.extend(norm_mu_error)
            total_max_error.extend(norm_max_error)
        print("Chromosomes:", self.chromosomes)
        print("Resolutions:", self.resolutions)
        print("Normalizations:", self.norms)
        print("Total error across all normalization vectors, all chromosomes, all resolutions ", total_error)
        print("Mean error across all normalization vectors, all chromosomes, all resolutions ", np.mean(total_mu_error))
        print("Max error across all normalization vectors, all chromosomes, all resolutions ", np.max(total_max_error))


def test_equality(file1: str, file2: str, norms: list):
    cf = CompareFiles(file1, file2, norms)
    cf.check_every_chromosome_pair()
    cf.check_all_normalizations()


if __name__ == '__main__':
    test_equality(sys.argv[1], sys.argv[2], sys.argv[3].split(","))
