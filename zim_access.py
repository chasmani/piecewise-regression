from zimscan import Reader

with Reader(open('gutenberg/zim/health.stackexchange.com_en_all_2020-05.zim', 'rb')) as reader:
    for record in reader:
        data = record.read()
        print(data)

