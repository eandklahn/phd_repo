from PyPDF2 import PdfFileWriter, PdfFileReader, PdfFileMerger
import argparse
import sys
import os

parser = argparse.ArgumentParser(description='Perform merging and splitting of pdf-files')

parser.add_argument('-t', '--type',
                    default='merge',
                    required=False,
                    dest='TYPE',
                    choices=['merge', 'split', 'rotate'])

parser.add_argument('-f', '--file',
                    default=None,
                    required=False,
                    dest='SPLITFILE')
                    
parser.add_argument('rest', nargs=argparse.REMAINDER)

def splitPDFintoPages(filename):
    """
    Splits the PDF given as filename into the constituent pages
    """
    
    inputpdf = PdfFileReader(open(filename, "rb"))
    
    for i in range(inputpdf.numPages):
        output = PdfFileWriter()
        output.addPage(inputpdf.getPage(i))
        with open("document-page%s.pdf" % i, "wb") as outputStream:
            output.write(outputStream)
            

def mergePDFs(inputList, outputname):
    """
    Merges the pdf-files given in inputlist into a single file.
    Inputlist should be a list of file paths, and outputname should be a text string
    """
    
    outputPDF = PdfFileMerger()

    for item in inputList:
        if not item.endswith('.pdf'):
            item += '.pdf'
        content = PdfFileReader(open(item, 'rb'))
        outputPDF.append(content)
    
    if not outputname.endswith('.pdf'):
        outputname += '.pdf'
    
    outputPDF.write(outputname)

def rotate_pdf(rest_list):
    
    file_in_name = rest_list[0]
    file_out_name = file_in_name.split('.')[0]+'_rot.pdf'
    
    infile = open(file_in_name, 'rb')
    outfile = open(file_out_name, 'wb')
    
    pdf_in = PdfFileReader(infile)    
    pdf_out = PdfFileWriter()
    for i in range(pdf_in.numPages):
        pdf_out.addPage(pdf_in.getPage(i).rotateCounterClockwise(90))
    
    pdf_out.write(outfile)
    
    infile.close()
    outfile.close()
    
if __name__ == '__main__':
    
    ARGS = parser.parse_args()
    
    if ARGS.TYPE=='merge':
        mergePDFs(ARGS.rest, 'merged')
    elif ARGS.TYPE=='split':
        if ARGS.SPLITFILE:
            splitPDFintoPages(ARGS.SPLITFILE)
        else:
            print('No file given to split')
    elif ARGS.TYPE=='rotate':
        rotate_pdf(ARGS.rest)