"""
AWS lambda function that handles bounce and complaint notifications from AWS
SES.
Marks bounced or compained about addresses as "bad" in the DB.
"""

import json
from nephele2.rds.db_utils import DBUtils


def lambda_handler(event, _):
    """
    Parses the lambda event and determines if action should be taken.
    """
    if 'Records' in event:
        for record in event['Records']:
            if 'Sns' in record and 'Message' in record['Sns']:
                message = json.loads(record['Sns']['Message'])
                if message['notificationType'] == "Bounce":
                    # Do this if email bounced
                    email_list = message['bounce']['bouncedRecipients']
                    for bad_addr in email_list:
                        DBUtils.set_bad_email(bad_addr['emailAddress'])
                elif message['notificationType'] == "Complaint":
                    # Do this if we receive a complaint
                    email_list = message['complaint']['complainedRecipients']
                    for bad_addr in email_list:
                        DBUtils.set_bad_email(bad_addr['emailAddress'])
