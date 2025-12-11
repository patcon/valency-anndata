from io import StringIO
import anndata as ad
import pandas as pd
from .utils import parse_polis_source
from polis_client import PolisClient

def load(source: str):
    """
    Accepts any of the following:
    - https://pol.is/report/<report_id>
    - https://pol.is/<conversation_id>
    - https://<host>/report/<report_id>
    - https://<host>/<conversation_id>
    - <conversation_id> (starts with digit; base_url=pol.is)
    - <report_id>       (starts with 'r'; base_url=pol.is)

    Returns:
        AnnData object
    """
    adata = ad.AnnData(X=None)

    convo_meta = parse_polis_source(source)
    client = PolisClient(base_url=convo_meta.base_url)
    # client = PolisClient(base_url=convo_meta.base_url, xid="foobar")

    if convo_meta.report_id:
        votes_csv_text = client.get_export_file(
            filename="votes.csv",
            report_id=convo_meta.report_id,
        )
        votes_csv_df = pd.read_csv(StringIO(votes_csv_text))
        votes_csv_df.sort_values(by="timestamp", ascending=True, inplace=True)

        adata.uns["votes_csv"] = votes_csv_df

        report = client.get_report(report_id=convo_meta.report_id)
        if report:
            convo_meta.conversation_id = report.conversation_id

    elif convo_meta.conversation_id:
        votes_list = client.get_all_votes_slow(conversation_id=convo_meta.conversation_id)
        votes_api_df = pd.DataFrame(votes_list)
        votes_api_df.sort_values(by="modified", ascending=True, inplace=True)

        adata.uns["votes_api"] = votes_api_df

    else:
        raise

    statements = client.get_comments(conversation_id=convo_meta.conversation_id)
    adata.uns["statements"] = pd.DataFrame(statements)

    # if convo_meta.conversation_id:
    #     xids = client.get_xids(conversation_id=convo_meta.conversation_id)
    #     adata.uns["xids"] = pd.DataFrame(xids)

    return adata
