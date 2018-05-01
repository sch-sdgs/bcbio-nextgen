"""Access StarlIMS functionality via sdgs webservice.
"""
from suds.client import Client
from suds.sudsobject import asdict

class StarLimsApi:
    def __init__(self):
        pass

    def _get_client(self):
        """
        contacts the webservice and gets a client with suds

        :return: client object
        """
        url = "http://10.182.155.37/StarLimsWeb.asmx?WSDL"
        client = Client(url)
        return client

    def run_details(self, run_id):
        """
        method to get the patients on a sequencer run

        :param run_id: hiseq folder name i.e. 150709_D00461_0040_AHTC7VADXX
        :return: list of dicts containing patients on the run
        """
        runName = {"runName": run_id}
        client = self._get_client()
        response = client.service.NGSByRun(**runName)
        result = recursive_asdict(response)
        return self.format_result(result)

    def all_ngs_test(self):
        """
        method to get all NGS patients

        :param run_id: hiseq folder name i.e. 150709_D00461_0040_AHTC7VADXX
        :return: list of dicts containing all patients
        """
        client = self._get_client()
        response = client.service.NGSAnalysis()
        result = recursive_asdict(response)
        return self.format_result(result)

    def ngs_by_container(self, container_id):
        """
        method to get patient data based on dna container id

        :param container_id: container id i.e. S1604013-02
        :return: patient info
        """
        runName = {"containerID": container_id}
        client = self._get_client()
        response = client.service.NGSByContainer(**runName)
        result = recursive_asdict(response)
        return self.format_result(result)[0]

    def format_result(self,result):
        """
        format the dict suds object and add patients to a list

        :param result:
        :return: list of patient dicts
        """
        out = []
        try:
            for i in result["diffgram"][0]["DocumentElement"][0]["results"]:
                info = {}
                for k, v in i.items():
                    info[k] = "".join(v)
                out.append(info)
            return out
        except TypeError:
            return []



def recursive_asdict(d):
    """
    Convert Suds object into serializable format.

    :param d: suds object
    :return: dict of suds object
    """
    out = {}
    for k, v in asdict(d).iteritems():
        if hasattr(v, '__keylist__'):
            out[k] = recursive_asdict(v)
        elif isinstance(v, list):
            out[k] = []
            for item in v:
                if hasattr(item, '__keylist__'):
                    out[k].append(recursive_asdict(item))
                else:
                    out[k].append(item)
        else:
            out[k] = v
    return out
