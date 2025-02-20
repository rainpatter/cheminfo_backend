import requests

BASE = "http://127.0.0.1:5000/"

# response = requests.get(BASE + 'helloworld')
response = requests.post(BASE + 'helloworld/CC=O')
print(response.json())