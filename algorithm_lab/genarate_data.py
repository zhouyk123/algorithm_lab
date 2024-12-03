import os
from python import mycalc

if __name__ == "__main__":
    with open('query.txt') as f:
        query = f.read().split('\n')[0]
    with open('ref.txt') as f:
        ref = f.read().split('\n')[0]
    os.system("g++ -O3 -o wqy.exe wqy.cpp")
    for i in range(-30,-3):
        for j in range(-30, -3):
            os.system("./wqy.exe %d %d" % (i, j))
            with open("result.txt") as f:
                result = f.read().split('\n')[0]
            result = result[1:-1].split(',')
            result = [int(x) for x in result]
            f.close()
            with open("args.txt", "a") as f:
                print(i, j, mycalc(result, ref, query, 0))
                f.write("%d %d %d\n" % (i, j, mycalc(result, ref, query, 0)))
            f.close()