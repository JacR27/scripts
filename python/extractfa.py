import sys

def main():
    seq = ""
    for i in sys.stdin:
        if not i.startswith(">"):
            seq += i.strip()
            
    seq = seq[154575326:154575327+160]
    ali = "TTTGGGGTGTGGGTTGGTGTGTGTGGTGTGTGTGTGTGTGTGTGGTGTGTGTGTGATGTGTGTGTATGTGTGTGATGTGTGTGTAGTGTGTGTGTGGTGTTTGTGTATGTGTGTGTGGTGTGTATGTGTGTGGTGTGTGTATGTGTGTGN"
    ali = ali[17:]
    for i,base in enumerate(ali):
        if base != seq[i]:
            print(str(i),base ,seq[i])
            
main()
            
