#!/usr/bin/env python3.6

import pandas as pd
import os
import argparse
from datetime import datetime, date


def main():

    # Parse argument
    parser = argparse.ArgumentParser(description="This software generates an html variants report.")
    parser.add_argument("-s", "--sampleid",help="Sample ID",dest="sampleid", type=str, required = True)
    parser.add_argument("-w", "--workdir", help="This is the base work directory.", dest="workdir", type=str, required = True)
    args = parser.parse_args()

    today = date.today()

    # sampleid = "BC05103"
    # workdir="../results/BC05103/confident_set/"
    print('[generate_report.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Generating variants report')

    files = os.listdir(workdir)
    files = [file for file in files if 'report.html' not in file]

    html_string = '''
    <html>
        <head>
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
            <style>body{ margin:0 100; background:whitesmoke; }</style>
        </head>
        <body>
        '''

    html_string = html_string + "<h1>" + sampleid + " Variants Report</h1>" + "<h2>Generated on: " + str(today) + "<h2/>"

    for idx, file in enumerate(files):
        try:
            df = pd.read_csv(workdir+file, sep='\t')
            summary_table = df.to_html().replace('<table border="1" class="dataframe">', '<table class="table table-striped">')
            html_string = html_string + "<h3>" + file + "</h3>"
            html_string = html_string + summary_table

        except:
            pass

    html_string = html_string + "</body>"
    html_string = html_string + "</html>"

    # Write the html report
    f = open(workdir + '/report.html','w')
    f.write(html_string)
    f.close()

    print('[generate_report.py]:  ' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ' Pipeline finished successfully')


if __name__=="__main__":
    main()
