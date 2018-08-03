pipeline {
    agent none
    stages {
        stage('Build') {
            agent {
                docker {
                    image 'quay.io/thoba/pythongit'
                }
            }
            steps {
                echo 'Building...'
                sh 'pip install -r requirements.txt'
                sh 'pip install -e .'
            }
        }
        stage('Test') {
            agent {
                docker {
                    image 'quay.io/thoba/pythongit'
                }
            }
            steps {
                echo 'Testing...'
                sh 'pip install -r requirements.txt'
                sh 'pytest -v test/ '
            }
        }
    }
}