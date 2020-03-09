from numpy import arange, array, hstack, int8, vstack, zeros


atomic_nums = vstack((hstack((array([1]), zeros(30), array([2]))),
                      hstack((arange(3, 5), zeros(24), arange(5, 11))),
                      hstack((arange(11, 13), zeros(24), arange(13, 19))),
                      hstack((arange(19, 22), zeros(14), arange(22, 37))),
                      hstack((arange(37, 40), zeros(14), arange(40, 55))),
                      arange(55, 87),
                      arange(87, 119))).astype(int8)

if __name__ == '__main__':
    from pandas import DataFrame
    print('Periodic table (atomic numbers, long form)')
    print(DataFrame(atomic_nums))

