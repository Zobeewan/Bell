"""Construit le PDF du papier depuis la source Markdown, avec equations rendues."""
from pathlib import Path
import hashlib
import html
import re
import json
import matplotlib
matplotlib.use('Agg')
from matplotlib import mathtext, font_manager
from PIL import Image as PILImage
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY,TA_LEFT,TA_CENTER
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Image,Table,TableStyle,PageBreak,KeepTogether,Preformatted

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'paper/Bell_Pair_Coherence_v2.md'
DEST=SOURCE.with_suffix('.pdf')
TMP=ROOT/'tmp/pdfs/bell_v2'; TMP.mkdir(parents=True,exist_ok=True)
FONTDIR=Path(matplotlib.get_data_path())/'fonts/ttf'
for name,file in [('Body','DejaVuSerif.ttf'),('BodyBold','DejaVuSerif-Bold.ttf'),('BodyItalic','DejaVuSerif-Italic.ttf'),
                  ('Sans','DejaVuSans.ttf'),('SansBold','DejaVuSans-Bold.ttf'),('Mono','DejaVuSansMono.ttf')]:
    pdfmetrics.registerFont(TTFont(name,str(FONTDIR/file)))
pdfmetrics.registerFontFamily('Body',normal='Body',bold='BodyBold',italic='BodyItalic',boldItalic='BodyBold')
pdfmetrics.registerFontFamily('Sans',normal='Sans',bold='SansBold',italic='Sans',boldItalic='SansBold')
WIDTH=A4[0]-104
styles={
 'body':ParagraphStyle('body',fontName='Body',fontSize=10.8,leading=14.8,spaceAfter=7,alignment=TA_JUSTIFY),
 'title':ParagraphStyle('title',fontName='SansBold',fontSize=23,leading=29,spaceAfter=14,textColor=colors.HexColor('#182c41')),
 'h2':ParagraphStyle('h2',fontName='SansBold',fontSize=14.2,leading=19,spaceBefore=8,spaceAfter=10,keepWithNext=True,textColor=colors.HexColor('#182c41')),
 'h3':ParagraphStyle('h3',fontName='SansBold',fontSize=11.2,leading=15,spaceBefore=6,spaceAfter=7,keepWithNext=True),
 'caption':ParagraphStyle('caption',fontName='BodyItalic',fontSize=9.3,leading=12.8,spaceAfter=15),
 'table':ParagraphStyle('table',fontName='Sans',fontSize=9.0,leading=12),
 'meta':ParagraphStyle('meta',fontName='Sans',fontSize=10,leading=14,spaceAfter=8,textColor=colors.HexColor('#5b6570')),
 'refs':ParagraphStyle('refs',fontName='Body',fontSize=9.5,leading=13.2,spaceAfter=9),
 'code':ParagraphStyle('code',fontName='Mono',fontSize=8.2,leading=12,spaceAfter=9),
}

def inline(text):
    text=text.replace('\u2013','-').replace('\u2014','-').replace('〉','⟩')
    text=html.escape(text)
    text=re.sub(r'\[([^\]]+)\]\((https?://[^\s)]+)\)',lambda m:f'<link href="{m[2]}" color="#205c86"><u>{m[1]}</u></link>',text)
    text=re.sub(r'\*\*([^*]+)\*\*',r'<b>\1</b>',text)
    return text

def equation(latex):
    latex=' '.join(latex.split())
    f=TMP/('eq_'+hashlib.sha256(latex.encode()).hexdigest()[:12]+'.png')
    if not f.exists():
        mathtext.math_to_image('$'+latex+'$',str(f),prop=font_manager.FontProperties(size=13),dpi=300,format='png')
    w,h=PILImage.open(f).size
    width=w*72/300; height=h*72/300
    scale=min(1.,WIDTH/width)
    im=Image(str(f),width=width*scale,height=height*scale)
    im.hAlign='CENTER'
    return KeepTogether([Spacer(1,5),im,Spacer(1,12)])

def parse_page(text,index):
    lines=text.strip().splitlines(); result=[]; i=0
    references='## Références' in text
    while i<len(lines):
        line=lines[i].strip()
        if not line: i+=1; continue
        if line=='$$':
            j=i+1
            while j<len(lines) and lines[j].strip()!='$$': j+=1
            result.append(equation('\n'.join(lines[i+1:j]))); i=j+1; continue
        if line.startswith('~~~'):
            j=i+1
            while j<len(lines) and not lines[j].startswith('~~~'): j+=1
            result.append(Preformatted('\n'.join(lines[i+1:j]),styles['code'])); i=j+1; continue
        if line.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].strip().startswith('|'):
                cells=[x.strip() for x in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch('[-: ]+',x) for x in cells):
                    rows.append([Paragraph(inline(x),styles['table']) for x in cells])
                i+=1
            n=len(rows[0]); widths=([WIDTH*.44,WIDTH*.56] if n==2 else [WIDTH/n]*n)
            table=Table(rows,colWidths=widths,repeatRows=1,hAlign='LEFT')
            table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#eef2f5')),
                ('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#8998a7')),
                ('LINEBELOW',(0,-1),(-1,-1),.5,colors.HexColor('#8998a7')),
                ('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),7),
                ('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6)]))
            result.extend([table,Spacer(1,12)]); continue
        if line.startswith('!['):
            m=re.match(r'!\[([^\]]*)\]\((.+)\)',line)
            path=(SOURCE.parent/m[2]).resolve(); w,h=PILImage.open(path).size
            im=Image(str(path),width=WIDTH,height=WIDTH*h/w)
            im.hAlign='CENTER'; result.append(im); result.append(Spacer(1,4)); i+=1; continue
        if line.startswith('# '): result.append(Paragraph(inline(line[2:]),styles['title'])); i+=1; continue
        if line.startswith('## '): result.append(Paragraph(inline(line[3:]),styles['h2'])); i+=1; continue
        if line.startswith('### '): result.append(Paragraph(inline(line[4:]),styles['h3'])); i+=1; continue
        if line.startswith('*') and line.endswith('*'):
            result.append(Paragraph(inline(line[1:-1]),styles['caption'])); i+=1; continue
        if re.match(r'^\d\. ',line):
            result.append(Paragraph(inline(line),styles['body'])); i+=1; continue
        para=[line]; i+=1
        while i<len(lines) and lines[i].strip() and not lines[i].startswith(('#','!','|','$$','~~~')):
            para.append(lines[i].strip()); i+=1
        text=' '.join(para)
        kind='refs' if references else 'body'
        if index==0 and (text.startswith('Version 2') or text.startswith('17 septembre')): kind='meta'
        result.append(Paragraph(inline(text),styles[kind]))
    return result

def footer(canvas,doc):
    canvas.saveState(); canvas.setFillColor(colors.HexColor('#6d7781')); canvas.setFont('Sans',8)
    canvas.drawString(52,27,'Cohérence de paire et corrélations de Bell | v2.1')
    canvas.drawRightString(A4[0]-52,27,str(doc.page))
    canvas.restoreState()

pages=SOURCE.read_text(encoding='utf-8').split('<!-- page -->')
story=[]
for i,page in enumerate(pages):
    if i: story.append(PageBreak())
    story.extend(parse_page(page,i))
doc=SimpleDocTemplate(str(DEST),pagesize=A4,leftMargin=52,rightMargin=52,topMargin=45,bottomMargin=46,
                      title='Cohérence de paire et corrélations de Bell',author='',
                      subject='Version 2.1 de Bell Indistinguishability')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
from pypdf import PdfReader
r=PdfReader(DEST)
report=dict(pages=len(r.pages),planned_sections=len(pages),links=sum(len(p.get('/Annots',[])) for p in r.pages),
            characters_by_page=[len(p.extract_text()) for p in r.pages])
(TMP/'build_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(str(DEST)); print(report)
