import sys
import re

def main():

    fn = 'pypeit_cite.bib'

    lines = read_bibtext(fn)

    years, months, refereed = parse_lines(lines)
    write_out(years, months, refereed)

def read_bibtext(fn):

    fp = open(fn)
    lines = fp.readlines()
    fp.close()

    return lines

def parse_tag(line):

    data = line.split()
    tag = data[0]

    values = data[2:]

    value = ''
    if len(values) > 1:
        value = " ".join(values)
    elif len(values) == 1:
        value = values[0]
    else:
        print(line)
    

    return tag, value

def get_month(value):

    months = {
        'jan' : 0,
        'feb' : 1,
        'mar' : 2,
        'apr' : 3,
        'may' : 4,
        'jun' : 5,
        'jul' : 6,
        'aug' : 7,
        'sep' : 8,
        'oct' : 9,
        'nov' : 10,
        'dec' : 11
        }

    return months[value]

def parse_lines(lines):

    in_entry = False

    years = []
    months = []
    refereed = []
    
    for line in lines:
        line = line.strip()
        if in_entry:
            if line[0] == "}":
                in_entry = False
            else:
                tag, value = parse_tag(line[0:-1])

                if tag == 'year':
                    years.append(value)
                elif tag == 'journal':
                    m = re.search("arXiv", value)
                    if m:
                        refereed.append(0)
                    else:
                        refereed.append(1)
                elif tag == 'month':
                    months.append(get_month(value))

        else:
            if len(line) > 0 and line[0] == "@":
                in_entry = True

    
    return years, months, refereed

def write_out(years, months, refereed):

    ofp=open("pypeit_sum.csv", "w+")

    ofp.write("years,months,refereed\n")
    
    for i in range(0, len(years)):

        oline = f"{years[i]},{months[i]},{refereed[i]}\n"
        ofp.write(oline)

    ofp.close()

    
if __name__ == "__main__":
    main()

    
