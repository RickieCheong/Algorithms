def z_algo(string):
    arr = list(string)
    l, r = 0, 0
    box = [0] * len(string)
    box[0] = len(string)
    length = len(string)

    for i in range(1, len(string)):
        # first case, checks to see if the index is beyond the z_box
        # Base case
        if i > r:
            l, r = i, i

            # ignores repetitive comparison and just assign value based on previously repeated values
            while r < length and arr[r - l] == arr[r]:
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
                # new pointer for the z box on the left side
                l = i
                while r + 1 < length and arr[r - l + 1] == arr[r + 1]:
                    r += 1

                box[i] = r - l + 1
    return box