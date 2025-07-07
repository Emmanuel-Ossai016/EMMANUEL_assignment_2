rule hello:
    output:
        "hello.txt"
    shell:
        "echo Hello World > hello.txt"