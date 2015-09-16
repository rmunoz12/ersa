import os
import math
from PIL import Image, ImageDraw, ImageFont
from visualization.chromsegments import *

def draw_rounded_rectangle(input, x1, y1, x2, y2, fillcolor = 'lightblue'):
    total_length = x2 - x1
    total_width = y2 - y1
    rounded_radius = total_width / 2

    draw = ImageDraw.Draw(input)
    draw.rectangle((x1, y1 + 2, x2, y2 - 2), fill = fillcolor)
    #draw.rectangle((x1 + rounded_radius, y1, x2 - rounded_radius, y2), fill = fillcolor)
    #draw.ellipse((x1, y1, x1 + 2 * rounded_radius, y2), fill = fillcolor)
    #draw.ellipse((x2 - 2 * rounded_radius, y1, x2, y2), fill = fillcolor)

def draw_chromosome(input, x1, y1, x2, y2, centromere_frac = 0.5, fillcolor = 'lightblue'):
    cent_x = centromere_frac * (x2 - x1)
    draw_rounded_rectangle(input, x1, y1, cent_x, y2, fillcolor)
    draw_rounded_rectangle(input, cent_x, y1, x2, y2, fillcolor)

def round_edges(input, x1, y1, x2, y2, centromere_frac = 0.5, fillcolor = 'white'):
    draw = ImageDraw.Draw(input)
    radius = (y2 - y1) / 2
    for i in range(0, int(y2 - y1)):
        line_len = radius - math.sqrt(max(2 * radius * i - i * i, 0))
        if line_len >= 1:
            draw.line((x1, y1 + i, x1 + line_len, y1 + i), fill = fillcolor)
            draw.line((centromere_frac * (x2 - x1) - line_len, y1 + i, centromere_frac * (x2 - x1) + line_len, y1 + i), fill = fillcolor)
            draw.line((x2 - line_len, y1 + i, x2, y1 + i), fill = fillcolor)

def create_pretty_chrom_plot(pair_dict, chrom_dict, xdim = 1024, ydim = 1024, pair = "unspecified"):
    # select pair if unspecified
    if pair == "unspecified":
        pair = random.choice(list(pair_dict.keys()))
    elif pair not in pair_dict:
        pair = pair.split(':')[1] + ':' + pair.split(':')[0]
        print(pair)
        if pair not in pair_dict:
            print("No matching entry")
            return

    # set up font
    fontfile = os.path.join(os.path.dirname(__file__), 'arial.ttf')
    print(fontfile)
    chrom_font = ImageFont.truetype(fontfile, 25)

    input = Image.new('RGB', (xdim, ydim), 'white')
    draw = ImageDraw.Draw(input)

    num_autosomes = 22
    chrom_sep = 10
    chrom_width = (ydim - 2 * chrom_sep) / 22 - chrom_sep
    chrom_longest_len = xdim * 0.8

    # get longest chromosome length by which to scale the drawn length
    longest_chrom = 0
    for chrom in range(1, num_autosomes + 1):
        if chrom_dict[chrom][0] > longest_chrom:
            longest_chrom = chrom_dict[chrom][0]

    print(pair)
    for chrom in range(1, num_autosomes + 1):
        x1 = 10
        x2 = x1 + chrom_dict[chrom][0] / longest_chrom * chrom_longest_len
        y1 = chrom_sep + (chrom_sep + chrom_width) * (chrom - 1)
        y2 = y1 + chrom_width
        cent_frac = chrom_dict[chrom][1] * 1e6 / chrom_dict[chrom][0]

        # draw chromosome
        draw_chromosome(input, x1, y1, x2, y2, cent_frac)

        # print chrom label
        draw.text((x2 + chrom_sep/2, y1), str(chrom), 'black', font = chrom_font)

        # draw shared segments on chromosomes
        segs = []
        for seg in [i for i in pair_dict[pair] if i[0] == chrom]:
            total_len = chrom_dict[chrom][0]
            print("chrom " + str(chrom) + " (" + str(total_len) + "):", seg[1], "-", seg[2], "|", seg[1]/total_len, "-", seg[2]/total_len)
            #segs.append((seg[1], seg[2] - seg[1])) #start, width
            seg_start = x1 + seg[1] / total_len * (x2 - x1)
            seg_end = x1 + seg[2] / total_len * (x2 - x1)


             # make this rounded at the edges when required
            draw_rounded_rectangle(input, seg_start, y1, seg_end, y2, 'blue')
        round_edges(input, x1, y1, x2, y2, cent_frac)
    input.show()

def pretty_parse_and_plot(path, t, h):
    pair_dict = get_pair_dict(path)
    chrom_dict = get_chrom_lengths('autosomelengths.txt')
    create_pretty_chrom_plot(pair_dict, chrom_dict)


def main():
    args = sys.argv

    h = 10  #high thres (cM)
    t = 2.5
    pretty_parse_and_plot(args[1], t, h)

if __name__ == '__main__':
    main()