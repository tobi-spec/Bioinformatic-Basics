from random import randint


def random_DNA(lenght):

    random_DNA_list =[]
    nucleotides = ["A", "T", "G", "C"]
    for i in range(0,lenght):
        random = randint(0,3)
        random_DNA_list.append(nucleotides[random])
    
    return "".join(random_DNA_list)

