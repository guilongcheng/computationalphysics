# -*- coding: utf-8 -*-


def TranName(names):
    new_names = ''
    for name in names:
        split_names = name.split()
        len_name = len(split_names)
        tmp = ''
        for i in range(len_name-1):
            tmp = tmp +split_names[i][0].upper() + '.~'
        tmp = tmp + split_names[-1].capitalize()
        new_names = new_names + tmp + ", "
    new_names = new_names[:-2]
    print("new name format as ",new_names)

def main():
    author = ["Gerg P. Engel", "C. B. Lang", "Daniel Mohler", "Andreas Schafer"]
    TranName(author)


if __name__ == "__main__":
    main()

