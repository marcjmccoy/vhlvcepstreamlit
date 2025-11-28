import json
from collections import defaultdict
import requests

API_KEY = "***"
API_URL = "https://api.perplexity.ai/chat/completions"
HEADERS = {
    "Authorization": f"Bearer {API_KEY}",
    "Content-Type": "application/json",
}

def query_llm_summarize(manager_node, subgraph_json):

    """
    Query Perplexity LLM to extract creative introduction paths based on local graph structure.
    """

    prompt = (
        f"Given this manager or organization node: {json.dumps(manager_node, ensure_ascii=False)}.\n"
        f"Relevant graph neighborhood (entities and links):\n"
        f"{json.dumps(subgraph_json, ensure_ascii=False)}\n"
        "***"
    )

    data = {
        "model": "sonar-pro",
        "messages": [
            {"role": "system", "content": "***"},
            {"role": "user", "content": prompt}
        ],
        "max_tokens": 500,
        "temperature": 0.35,
    }

    response = requests.post(API_URL, headers=HEADERS, json=data)
    response.raise_for_status()
    result = response.json()
    answer = result.get("choices", [{}])[0].get("message", {}).get("content", "").strip()
    return answer