"""Build the README examples from private TIFFs and published coordinates."""
import argparse
import hashlib
from pathlib import Path
import textwrap
import numpy as np
import tifffile
from PIL import Image, ImageDraw, ImageFont
from analyze import read_csv, detections
from atlas import ontology, display_group
from config import ROOT, RESULTS

BLUE = '#00bfff'
PINK = '#ff66cc'
GOLD = '#ffdf43'


def font(size):
    for name in ['DejaVuSans.ttf', '/System/Library/Fonts/Helvetica.ttc']:
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            pass
    return ImageFont.load_default(size=size)


def background(raw, low, high):
    lut = np.uint8(np.clip((np.arange(65536, dtype=float)-low)/max(high-low, 1), 0, 1)**.65*255)
    return Image.fromarray(lut[raw]).convert('RGB')


def overlay(base, points, scale, color):
    panel = base.copy()
    draw = ImageDraw.Draw(panel)
    radius = max(1.3, 3.5*scale)
    for y, x in points:
        x, y = x*scale, y*scale
        draw.ellipse((x-radius, y-radius, x+radius, y+radius), outline=color, width=max(1, round(scale)))
    return panel


def whole_slice(raw, sec, low, high, output):
    a, _ = detections(sec)
    b, _ = detections(sec, 'alternative')
    width = 1000
    scale = width/raw.shape[1]
    height = round(raw.shape[0]*scale)
    base = background(raw, low, high).resize((width, height), Image.Resampling.LANCZOS)
    page = Image.new('RGB', (3000, height+160), 'white')
    draw = ImageDraw.Draw(page)
    draw.text((15, 10), f'Mouse 1 | section {sec} | whole-slice GFP detections', font=font(32), fill='black')
    for col, (label, points, color) in enumerate([('Raw', [], BLUE), (f'Model: {len(a):,}', a, BLUE), (f'Alternative: {len(b):,}', b, PINK)]):
        draw.text((col*width+15, 57), label, font=font(26), fill='black')
        page.paste(overlay(base, points, scale, color), (col*width, 100))
    draw.text((15, height+116), 'Both methods use the same atlas + tissue area. Dots mark detections; the alternative is not ground truth.', font=font(22), fill='black')
    page.save(output, optimize=True)


def counting_crop(raw, sec, um, row):
    y, x, size = int(row['y']), int(row['x']), int(row['size'])
    crop = raw[y:y+size, x:x+size]
    a, _ = detections(sec)
    b, _ = detections(sec, 'alternative')
    def local(points):
        return points[np.all((points >= [y, x]) & (points < [y+size, x+size]), axis=1)]-[y, x]
    a, b = local(a), local(b)
    assert len(a) == int(row['full']) and len(b) == int(row['alternative'])
    lo, hi = np.percentile(crop, [1, 99.7])
    base = background(crop, lo, hi).resize((640, 640), Image.Resampling.NEAREST)
    page = Image.new('RGB', (1980, 830), 'white')
    draw = ImageDraw.Draw(page)
    draw.text((15, 10), f'Section {sec} | local comparison | y={y}, x={x} | {size} x {size} pixels', font=font(28), fill='black')
    for col, (label, points, color) in enumerate([('Raw', [], BLUE), (f'Model: {len(a)}', a, BLUE), (f'Alternative: {len(b)}', b, PINK)]):
        draw.text((col*660+15, 58), label, font=font(26), fill='black')
        panel = overlay(base, points, 640/size, color)
        ink = ImageDraw.Draw(panel)
        ink.line((15, 615, 15+100/um*640/size, 615), fill='white', width=4)
        ink.text((15, 581), '100 um', font=font(22), fill='white', stroke_width=1, stroke_fill='black')
        page.paste(panel, (col*660+10, 105))
    draw.text((15, 763), 'Targeted disagreement example; identical contrast across methods. This is not an unbiased accuracy sample.', font=font(22), fill='black')
    return page


def atlas_figure(raw, labels, sec, um, low, high, output):
    nodes, paths = ontology()
    palette = {int(r['structure_id']): r['display_hex'] for r in read_csv(ROOT/'atlas/display_palette.csv')}
    hotspots = [r for r in read_csv(RESULTS/'hotspot_summary.csv') if int(r['section']) == sec]
    width = 1200
    scale = width/raw.shape[1]
    height = round(raw.shape[0]*scale)
    ys = np.minimum(((np.arange(height)+.5)/scale).astype(int), raw.shape[0]-1)
    xs = np.minimum(((np.arange(width)+.5)/scale).astype(int), raw.shape[1]-1)
    sample = labels[np.ix_(ys, xs)]
    group = np.zeros(sample.shape, dtype=np.int64)
    rgb = np.zeros((*sample.shape, 3), dtype=np.uint8)
    for sid in np.unique(sample):
        if sid == 0:
            continue
        gid = display_group(int(sid), nodes, paths)
        mask = sample == sid
        group[mask] = gid
        rgb[mask] = tuple(bytes.fromhex(palette[gid]))
    edge = np.zeros(sample.shape, dtype=bool)
    edge[1:] |= group[1:] != group[:-1]
    edge[:, 1:] |= group[:, 1:] != group[:, :-1]
    rgb[edge] = [220, 220, 220]
    atlas = Image.fromarray(rgb)
    base = background(raw, low, high).resize((width, height), Image.Resampling.LANCZOS)
    page = Image.new('RGB', (2440, height+305+120*((len(hotspots)+1)//2)), 'white')
    draw = ImageDraw.Draw(page)
    draw.text((20, 10), f'Mouse 1 | section {sec} | Allen locations of GFP hotspots', font=font(32), fill='black')
    draw.text((20, 55), 'Existing atlas alignment | Distinct region colors | Circles are exploratory review windows', font=font(25), fill='black')
    for col, (panel, title) in enumerate([(base, 'Raw slice'), (atlas, 'Aligned atlas; cortical layers pooled for display')]):
        ink = ImageDraw.Draw(panel)
        for row in hotspots:
            x, y = float(row['center_x_px'])*scale, float(row['center_y_px'])*scale
            rad = float(row['radius_um'])/um*scale
            if col == 1:
                ink.ellipse((x-rad-1, y-rad-1, x+rad+1, y+rad+1), outline='black', width=5)
            ink.ellipse((x-rad, y-rad, x+rad, y+rad), outline=GOLD, width=3)
            ink.text((max(2, x-rad), max(2, y-rad-23)), row['hotspot'], font=font(25), fill=GOLD, stroke_width=2, stroke_fill='black')
        draw.text((20+col*1220, 100), title, font=font(25), fill='black')
        page.paste(panel, (20+col*1220, 140))
    for k, row in enumerate(hotspots):
        x, y = 20+(k % 2)*1220, 160+height+(k//2)*120
        draw.rectangle((x, y+3, x+18, y+21), fill='#'+palette[int(row['dominant_area_id'])], outline='black')
        draw.text((x+28, y), f"{row['hotspot']}. {row['dominant_acronym']} — {row['dominant_name']}", font=font(25), fill='black')
        draw.text((x+28, y+31), f"{row['detections']} detections | dominant area {row['dominant_detections']}/{row['detections']} ({100*float(row['dominant_fraction']):.0f}%)", font=font(22), fill='black')
        for line, text in enumerate(textwrap.wrap(row['composition'], width=100)[:2]):
            draw.text((x+28, y+60+24*line), text, font=font(22), fill='#333333')
    draw.text((20, page.height-105), 'Mixed circles retain all regional assignments. Counts are model detections, not confirmed cells.', font=font(24), fill='black')
    draw.text((20, page.height-65), 'Atlas alignment has not been independently validated. Custom colors encode region identity only.', font=font(22), fill='black')
    page.save(output, optimize=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=ROOT/'data/slices.csv')
    parser.add_argument('--image-root', type=Path, default=ROOT/'data')
    parser.add_argument('--output', type=Path, default=ROOT/'figures')
    args = parser.parse_args()
    records = {int(r['section']): r for r in read_csv(args.manifest)}
    review = read_csv(RESULTS/'review_windows.csv')
    args.output.mkdir(parents=True, exist_ok=True)
    crops = []
    for sec in [12, 16]:
        rec = records[sec]
        rawpath, labelpath = args.image_root/rec['raw_file'], args.image_root/rec['atlas_file']
        for path, expected in [(rawpath, rec['raw_sha256']), (labelpath, rec['atlas_sha256'])]:
            if not path.exists():
                raise SystemExit(f'Private source TIFF required: {path}. See README for the manifest format.')
            assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, path
        raw, labels = tifffile.imread(rawpath), tifffile.imread(labelpath)
        assert raw.dtype == np.uint16 and raw.shape == labels.shape
        lo, hi, um = float(rec['display_min']), float(rec['display_max']), float(rec['um_per_pixel'])
        whole_slice(raw, sec, lo, hi, args.output/f'whole_slice_{sec:02d}.png')
        atlas_figure(raw, labels, sec, um, lo, hi, args.output/f'atlas_{sec:02d}.png')
        row = next(r for r in review if int(r['section']) == sec)
        crops.append(counting_crop(raw, sec, um, row))
    gallery = Image.new('RGB', (1980, sum(im.height for im in crops)), 'white')
    y = 0
    for panel in crops:
        gallery.paste(panel, (0, y))
        y += panel.height
    gallery.save(args.output/'counting_examples.png', optimize=True)
    print('Wrote five README figures; original TIFFs were read without modification.')
