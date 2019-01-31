#!/usr/bin/python

import sys

def main():
	if not sys.argv[1] and sys.argv[2]:
		sys.exit(1)
		print("Usage: readsFromContigs.py <k> <string>")

	k=int(sys.argv[1])
	s=str(sys.argv[2])

	keep_going=True
	start=0
	end=start+k

	print(len(s))
	print(end)
	while end <= len(s):
		print(s[start:end])
		start+=1
		end+=1


#Call main function
if __name__ == '__main__':
    main()
