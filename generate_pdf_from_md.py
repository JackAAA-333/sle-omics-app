from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
import textwrap
import os

MD = 'outputs/report.md'
OUT = 'outputs/report.pdf'
IMG_DIRS = ['outputs/report_figs','outputs_advanced']

def draw_text(c, text, x, y, max_width):
    lines = []
    for para in text.split('\n'):
        wrapped = textwrap.wrap(para, width=100)
        if not wrapped:
            lines.append('')
        else:
            lines.extend(wrapped)
    for line in lines:
        c.drawString(x, y, line)
        y -= 12
    return y

def main():
    c = canvas.Canvas(OUT, pagesize=A4)
    width, height = A4
    # title
    c.setFont('Helvetica-Bold', 16)
    c.drawCentredString(width/2, height-50, 'SLE 多组学整合报告')
    y = height-80
    c.setFont('Helvetica', 10)
    if os.path.exists(MD):
        with open(MD, 'r', encoding='utf-8') as f:
            md = f.read()
    else:
        md = '报告文件未找到.'
    y = draw_text(c, md[:8000], 40, y, width-80)
    c.showPage()
    # add images from fig dirs
    for d in IMG_DIRS:
        if not os.path.isdir(d):
            continue
        for fname in os.listdir(d):
            if fname.lower().endswith(('.png','.jpg','.jpeg')):
                path = os.path.join(d, fname)
                try:
                    img = ImageReader(path)
                    iw, ih = img.getSize()
                    scale = min((width-80)/iw, (height-120)/ih, 1)
                    w = iw*scale; h = ih*scale
                    c.drawImage(img, 40, height-100-h, width=w, height=h)
                    c.showPage()
                except Exception:
                    continue
    c.save()
    print('PDF written:', OUT)

if __name__=='__main__':
    main()
