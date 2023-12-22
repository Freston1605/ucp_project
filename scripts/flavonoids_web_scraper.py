import requests
from bs4 import BeautifulSoup
import csv
from os import path

def request_chebi_id(url, output="chebi_ids.csv", ignore=None, timeout=10):
    """Make a request to the specified URL, extract ChEBI IDs from the web page, and write them to a CSV file.

    Args:
        url (str): The URL to make the request to.
        output (str, optional): The path to the output CSV file. Defaults to "chebi_ids.csv".
        ignore (str, optional): The ChEBI ID to ignore. Defaults to None.
        timeout (int, optional): The maximum number of seconds to wait for the request. Defaults to 10.
    """
    # Make a request to the URL with a specified timeout
    try:
        html_response = requests.get(url, timeout=timeout)
    except requests.exceptions.RequestException as e:
        print(f"Failed to retrieve the page. Error: {e}")
        return

    # Check if the request was successful (status code 200)
    if html_response.status_code == 200:
        # Parse the HTML content with BeautifulSoup
        soup = BeautifulSoup(html_response.text, "html5lib")

        # Find all <a> tags containing ChEBI IDs
        chebi_links = soup.find_all("a", href=lambda x: x and "chebiId" in x)

        # Extract and filter ChEBI IDs
        chebi_id_list = [
            link.text.split(":")[1]
            for link in chebi_links
            if "CHEBI:" in link.text and link.text.split(":")[1] != ignore
        ]

        # Write the list of ChEBI IDs to a CSV file
        with open(output, "w", newline="", encoding='utf-8') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(["ChEBI_ID"])  # Write header
            csv_writer.writerows([[chebi_id] for chebi_id in chebi_id_list])

        print(f"ChEBI IDs written to: {output}")
    else:
        print(f"Failed to retrieve the page. Status code: {html_response.status_code}")

# Execution

# Intended to extract all flavonoids in the ChEBI database's web page.
# Ignores the ID for flavonoids (72544) since it repeats on every molecule that is classified as a flavonoid.

# Variables
OUTPUT_PATH = "data/output/csv"
OUTPUT_FILE = "chebi_ids.csv"
OUTPUT_FILE = path.join(OUTPUT_PATH, OUTPUT_FILE)
URL = "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=72544"
IGNORE = "72544"

request_chebi_id(URL, output=OUTPUT_FILE, ignore=IGNORE)
