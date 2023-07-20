import sys
'''
Author: Rickie Cheong Jun Weng
Email: rche0041@student.monash.edu
Student ID: 30897580
'''


'''
Helper functions for pre-processing for Boyer Moore algorithm
'''

# bad character table that is modified to take in wild cards
def bad_character(pattern):

    length = len(pattern)

    # used for indexing the ASCII values of the alphabets
    characters = 97

    # Bad character table initialized with 0 for all 26 alphabets
    BC_table = [[0] * length for _ in range(26)]

    # filling in the bad character table with the appropriate values
    for i in range(len(pattern) - 1, -1, -1):

        # Filling the bad character with the wildcard in every row and column after the original wildcard placement
        if (ord(pattern[i])) == 46:
            for j in range(len(BC_table)):
                BC_table[j][i] = i + 1
                index = i
                while index < len(pattern) - 1 and BC_table[j][index + 1] == 0:
                    BC_table[j][index + 1] = BC_table[j][index]
                    index += 1
        # For all other 26 alphabets
        else:
            BC_table[ord(pattern[i]) - characters][i] = i + 1
        index = i

        # filling the rest of the row
        while index < len(pattern) - 1 and ord(pattern[i]) != 46 and BC_table[ord(pattern[i]) - characters][index + 1] == 0:
            BC_table[ord(pattern[i]) - characters][index + 1] = BC_table[ord(pattern[i]) - characters][index]
            index += 1

    return BC_table


# Normal matched prefix algorithm implemented based on Lectures and tutorials
def matched_prefix(string):
    length = len(string)
    box = z_algo(string)
    box.append(0)
    for i in range(length - 1, 0, -1):
        if box[i] + i == length:
            continue
        elif box[i] + i != length and i != length - 1:
            box[i] = box[i + 1]
    return box


# Normal Good suffix implemented based on lectures and tutorials
def good_suffix(pattern):
    length_pat = len(pattern)
    box = reversed_Z_algo(pattern)
    goodsuffix = [0] * (length_pat + 1)
    for i in range(length_pat - 1):
        j = length_pat - box[i] + 1
        goodsuffix[j - 1] = i + 1
    return goodsuffix


# Improvised Z algorithm to accept wildcard as any comparisons to be true
def z_algo(string):
    l, r = 0, 0
    box = [0] * len(string)
    box[0] = len(string)
    length = len(string)

    for i in range(1, len(string)):
        # Base case, checks to see if the index is beyond the z_box
        # case 1a
        if i > r:
            l, r = i, i

            # ignores repetitive comparison and just assign value based on previously repeated values
            while r < length and (string[r - l] == string[r] or ord(string[r]) == 46 or ord(string[r - l]) == 46):
                r += 1
            box[i] = r - l
            r -= 1

        else:
            # checks to see k index that corresponds to the current value in the current zbox
            k = i - l

            # if the k index is lesser than the remaining amount of boxes in the zbox, assign the value in
            # k index to the current i index
            # case 2a
            if box[k] < r - i + 1:
                box[i] = box[k]

            # if the value at the k index is greater than the remainder, we assign the number of remaining boxes
            # in the i index value
            # case 2b
            elif box[k] > r - i + 1:
                box[i] = r - i + 1

            # unsure if the value past the zbox is the same, so we check after the zbox until we find non
            # matching values
            # case 2c
            elif box[k] == r - i + 1:
                # new pointer for the new z box to check for more comparisons
                l = i
                while r + 1 < length and string[r - l + 1] == string[r + 1]:
                    r += 1

                box[i] = r - l + 1
                r -= 1
    return box


# Reversed version of the normal Z algorithm to calculate z score from suffix that takes in wild card
def reversed_Z_algo(string):
    # Initializes the z values
    l, r = 0, 0
    box = [0] * len(string)
    box[-1] = len(string)
    length = -len(string)
    for i in range(2, len(string) + 1):
        i = i * -1

        # case 1
        if i < l:
            l, r = i, i
            while l > length - 1 and (string[l - r - 1] == string[l] or ord(string[l]) == 46 or ord(string[l - r - 1]) == 46):
                l -= 1

            box[i] = -(l-r)
            l += 1

        else:
            k = i - r - 1

            # case 2a
            if box[k] < -(l - i - 1):
                box[i] = box[k]

            # case 2b
            elif box[k] > -(l - i - 1):
                box[i] = -(l - i - 1)

            # case 2c
            elif box[k] == -(l - i - 1):
                r = i

                # new zbox for additional comparisons to check for matches
                while l - 1 >= length and string[l - r - 2] == string[l - 1]:
                    l -= 1
                box[i] = -(l - r - 1)
                l += 1
    return box

'''
Main Function

Improved Boyer Moore algorithm to take accept wild card
'''
def Improved_boyer_moore(text, pattern):
    f = open("output_q2.txt", "w")

    # Initializes bad character table, Good suffix, and matched prefix
    BC_table = bad_character(pattern)
    GS_table = good_suffix(pattern)
    MP_table = matched_prefix(pattern)

    coords = 0
    temp = 0
    BC = 0
    GSMP = 0
    value = 0
    stop = 0
    resume = 0
    lst =[]
    characters = 97
    length_text = len(text)
    length_pat = len(pattern)

    while coords < length_text:
        counter = 0

        # To ensure we do not go past the possible length of the text
        if len(pattern) > len(text):
            break

        #
        for i in range(length_pat - 1, -1, -1):
            if (coords + i) >= length_text:
                coords += i
                break

            # Galil's Optimization, skips past values that we are sure is a match
            if stop < coords + i < resume:
                counter += 1
                continue

            # if mismatch happens, bad character rule and good suffix are calculated and sees which one provides a
            # greater jump for skipping comparisons, not including wildcard as mismatch
            if pattern[i] != text[coords + i] and (ord(pattern[i]) != 46):
                BC = i - BC_table[ord(text[coords + i]) - characters][i] + 1

                # If good suffix is used, we know which characters will have a match, and save values to skip comparing
                if GS_table[i + 1] > 0:
                    counter = 0
                    GSMP = length_pat - GS_table[i + 1]
                    suffix_length = length_pat - i
                    stop = coords + i
                    resume = coords + i + suffix_length - 1

                # Matched prefix if good suffix value is 0
                else:
                    GSMP = MP_table[i + 1]
                    stop = -1
                    resume = -1
                value = max(BC, GSMP)

                # If all 3 tables provide 0, there is no possible jump and have to skip through everything
                if value == 0:
                    value = len(pattern)
                coords += value
                break

            # full match with wildcard, add the values to file
            elif i == 0:
                lst.append(coords + 1)
                f.write(f"{coords + 1}\n")
                coords += length_pat - MP_table[1]

        # resets the values every time new comparison is used
        counter += 1
        if counter >= 2:
            stop = -1
            resume = -1
            counter = 0
    f.close()
    return


def read_file(filename):
    mystring = ""
    myfile = open(filename, 'r')
    for line in myfile:
        mystring += line.strip()
    myfile.close()
    return mystring


if __name__ == '__main__':
    _, filename1, filename2 = sys.argv
    command = sys.argv[0]
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]

    Improved_boyer_moore(read_file(filename1), read_file(filename2))


