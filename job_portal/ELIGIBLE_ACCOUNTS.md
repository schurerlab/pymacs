# Portal eligibility roster

The live roster is `/home/jxs794/lab-job-portal/jobs/eligible_accounts.json` on Kyle. It is the compatibility document used before the portal checks SSH or offers SSH setup.

Each entry requires an exact first name and University of Miami directory email, for example `jxs794@miami.edu`. Do not use vanity or departmental addresses such as `jmschulz@med.miami.edu`.

To add a member, append one object to the `accounts` array:

```json
{"name": "FirstName", "email": "universityid@miami.edu"}
```

The roster grants only the right to claim/setup a portal profile. The member still needs working passwordless SSH and administrator approval with a portal password before they can sign in or submit jobs.
