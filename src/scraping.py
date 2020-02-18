r"""
scraping

Module containing web scraping functions.

"""
from bs4 import BeautifulSoup
from os import mkdir
from re import search as matched
from requests import get as get_url
from wget import download


def download_link_locations(url, out=None, regex=None):
    r"""
    Download link locations from specified URL.

    Parameters
    ----------
    url : str
        Input url.
    out : str, optional
        Target directory.
    regex : str, optional
        Regex filter.

    Returns
    -------
    None

    """
    html = get_url(url).text
    soup = BeautifulSoup(html, 'html.parser')
    for a in soup.findAll('a'):
        href = a.get('href')
        if matched(regex, href):
            download(url+href, out)


if __name__ == '__main__':
    print('IMPORT-ONLY MODULE', __doc__, end='', sep='\n')

