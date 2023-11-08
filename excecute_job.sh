#!/bin/bash
python KIDS_test.py
echo "Job has finished. Mailing now.."
python /data/typhon2/email_script.py --rea "rhysalfshaw@gmail.com" --message "Job has finised."
echo "Job has finished. Mailed."vi 