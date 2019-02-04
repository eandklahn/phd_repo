from PyPDF2 import PdfFileWriter, PdfFileReader, PdfFileMerger

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
            
if __name__ == '__main__':
    
    list = ['0','1','2','3','4','5']
           
    mergePDFs(list, 'rental')
    
    #splitPDFintoPages('5pages.pdf')